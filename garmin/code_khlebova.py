from typing import Optional, Tuple
from dataclasses import dataclass
from enum import IntEnum
import pandas as pd
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

# Extract RRI, speed and timestamp values
df = pd.read_csv('rri.csv.gz', compression='gzip')
rri_values = df.iloc[:, 0].values
rri_values = np.array(rri_values)
speed_df = pd.read_csv('speed.csv.gz', compression='gzip')
speed_timestamps = pd.to_datetime(speed_df['timestamp'])
speed_values = speed_df['speed'].values


# Create a lookup dict from seconds-precision timestamp to speed value
speed_lookup = {}
for ts, speed in zip(speed_timestamps, speed_values):
    ts_seconds = ts.floor('s')
    speed_lookup[ts_seconds] = speed

# Create timestamps based on RRI values (cumulative time from start)
# RRI values are in milliseconds, cumulative sum gives time of each measurement
start_time = datetime(2025, 12, 1, 9, 0, 0)
cumulative_ms = np.cumsum(np.concatenate([[0], rri_values[:-1]]))
timestamps = [start_time + timedelta(milliseconds=float(ms))
              for ms in cumulative_ms]
timestamps = np.array(timestamps)

"""
Clean RRI outliers by replacing them with the median value.
RRIs below 300 ms (tachycardia or myoelectricity artifacts) and
above 3000 ms (sinus arrest or no signal) are considered outliers.
"""
rri_values[(rri_values >= 3000) | (rri_values <= 300)] = np.median(rri_values)

# RMSSD and HR are calculated over 5 min intervals


def calculate_rmssd(rri_window):
    if len(rri_window) < 2:
        return 0
    diff = np.diff(rri_window)
    return np.sqrt(np.mean(diff ** 2))


def calcualte_hr(rri_window):
    if len(rri_window) < 2:
        return 0
    return 60000 / np.mean(rri_window)


window_duration = timedelta(minutes=5)
rmssd_values = []
rmssd_timestamps = []
hr_values = []
hr_timestamps = []

for i in range(len(rri_values)):
    current_time = timestamps[i]
    window_start_time = current_time - window_duration

    # Find all RRI values within the 5-minute window
    # Get indices where timestamp is within [window_start_time, current_time]
    mask = (timestamps >= window_start_time) & (timestamps <= current_time)
    window = rri_values[mask]

    if len(window) >= 2:
        rmssd_values.append(calculate_rmssd(window))
        rmssd_timestamps.append(current_time)
        hr_values.append(calcualte_hr(window))
        hr_timestamps.append(current_time)
    else:
        rmssd_values.append(0)
        rmssd_timestamps.append(current_time)
        hr_values.append(0)
        hr_timestamps.append(current_time)

rmssd_values = np.array(rmssd_values)
rmssd_timestamps = np.array(rmssd_timestamps)
hr_values = np.array(hr_values)
hr_timestamps = np.array(hr_timestamps)

# Find the lactate threshold


class Status(IntEnum):
    """Return status for update_lactate_threshold function"""
    OK = 0                    # Processing continues normally
    SKIP_ZERO_SPEED = 1       # Skipped due to zero speed
    THRESHOLD_FOUND = 2       # Lactate threshold point found
    NO_CHANGE = 3             # No significant change detected


@dataclass
class LactateThresholdState:
    """State structure for lactate threshold detection"""
    # Trend detection parameters
    min_trend_length: int = 10

    # Previous values
    prev_rmssd: Optional[float] = None
    prev_index: int = 0

    # Decreasing trend tracking
    decrease_trend_count: int = 0
    decrease_trend_started: bool = False
    decrease_trend_confirmed: bool = False

    # Increasing trend tracking
    increase_trend_count: int = 0
    increase_trend_started: bool = False
    increase_start_index: Optional[int] = None
    increase_start_rmssd: Optional[float] = None
    increase_start_speed: Optional[float] = None


def update_lactate_threshold(
    rmssd: float,
    running_speed: float,
    heart_rate: int,
    state: LactateThresholdState
) -> Tuple[Status, LactateThresholdState]:
    """
    Update lactate threshold detection state with a new data point.

    Args:
        rmssd: Current RMSSD value in ms
        running_speed: Current running speed
        heart_rate: Current heart rate
        state: Current state of the detection algorithm

    Returns:
        Tuple of (status, updated_state)
    """
    # Skip if speed is zero
    if running_speed == 0:
        state.prev_index += 1
        return Status.SKIP_ZERO_SPEED, state

    # Need at least one previous value to compare
    if state.prev_rmssd is None:
        state.prev_rmssd = rmssd
        state.prev_index += 1
        return Status.OK, state

    current_val = round(rmssd, ndigits=2)
    prev_val = round(state.prev_rmssd, ndigits=2)
    current_index = state.prev_index

    # Phase 1: Looking for decreasing trend
    if not state.decrease_trend_confirmed:
        if prev_val > current_val:
            if not state.decrease_trend_started:
                state.decrease_trend_started = True
            state.decrease_trend_count += 1
        elif prev_val < current_val:
            if state.decrease_trend_started and state.decrease_trend_count >= state.min_trend_length:
                # Decreasing trend confirmed, now start looking for increasing trend
                state.decrease_trend_confirmed = True
                state.increase_trend_started = True
                state.increase_start_index = current_index - 1
                state.increase_start_rmssd = state.prev_rmssd
                state.increase_trend_count = 1
            else:
                # Reset - trend wasn't long enough
                state.decrease_trend_count = 0
                state.decrease_trend_started = False

    # Phase 2: Looking for increasing trend after decreasing trend
    else:
        if prev_val < current_val:
            if not state.increase_trend_started:
                state.increase_trend_started = True
                state.increase_start_speed = running_speed
                state.increase_start_index = current_index - 1
                state.increase_start_rmssd = state.prev_rmssd
            state.increase_trend_count += 1
        elif prev_val > current_val:
            if state.increase_trend_started and state.increase_trend_count >= state.min_trend_length:
                state.prev_rmssd = rmssd
                state.prev_index += 1
                return Status.THRESHOLD_FOUND, state
            else:
                # Reset increasing trend - wasn't long enough
                state.increase_trend_count = 0
                state.increase_trend_started = False
                state.increase_start_index = None
                state.increase_start_rmssd = None
                state.increase_start_speed = None

    state.prev_rmssd = rmssd
    state.prev_index += 1

    return Status.OK, state


state = LactateThresholdState(min_trend_length=10)

# Process all data points
for i in range(len(rmssd_values)):
    ts_seconds = pd.Timestamp(rmssd_timestamps[i]).floor('s')
    speed = speed_lookup.get(ts_seconds, 0.0)
    if np.isnan(speed):
        speed = 0.0

    status, state = update_lactate_threshold(
        rmssd=rmssd_values[i],
        running_speed=speed,
        heart_rate=hr_values[i],
        state=state
    )

    if status == Status.THRESHOLD_FOUND:
        print(f"Lactate threshold found at the following point:")
        print(f"  RMSSD: {state.increase_start_rmssd:.2f} ms")
        print(f"  Speed: {state.increase_start_speed}")
        break

# Create the plot with dual y-axes (RMSSD on left, speed on right)
fig, ax1 = plt.subplots(figsize=(12, 6))
# Plot RMSSD values against time on left y-axis
ax1.plot(rmssd_timestamps, rmssd_values, 'g-',
         linewidth=1.5, label='RMSSD values')
ax1.set_xlabel('Time')
ax1.set_ylabel('RMSSD (ms)', color='g')
ax1.tick_params(axis='y', labelcolor='g')
ax1.grid(True, alpha=0.3)

# Plot speed values on right y-axis
ax2 = ax1.twinx()
ax2.plot(speed_timestamps, np.array(speed_values),
         'purple', linewidth=1, label='Speed')
ax2.set_ylabel('Speed', color='purple')
ax2.tick_params(axis='y', labelcolor='purple')


# Add vertical and horizontal lines for the found lactate threshold
if state.increase_start_index is not None:
    found_time = rmssd_timestamps[state.increase_start_index]
    print(f"Found time: {found_time}")
    ts_seconds = pd.Timestamp(found_time).floor('s')
    # Vertical line for time at the lactate threshold
    ax1.axvline(x=found_time, color='r', linestyle='--',
                linewidth=1.5, label=f'Lactate threshold')
    # Horizontal line for RMSSD at the lactate threshold
    ax1.axhline(y=state.increase_start_rmssd, color='orange',
                linewidth=1.5, label=f'RMSSD at LT: {state.increase_start_rmssd:.1f}')

    ax1.annotate(f'LT time {found_time.hour}:{found_time.minute} and corresponding speed: {state.increase_start_speed:.2f} m/s', xy=(found_time, state.increase_start_rmssd),
                 xytext=(10, 10), textcoords='offset points', color='r')

ax1.legend(loc='upper left')
ax2.legend(loc=0)

plt.title('RMSSD Values over Time (with Speed on top axis)')
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
fig.autofmt_xdate()
plt.tight_layout()
plt.savefig('rri_plot.png', dpi=150)
plt.show()

# Create Plot 2: RRI, HR, RMSSD, and Speed timeseries
fig2, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

axes[0].plot(timestamps, rri_values, 'blue', linewidth=1)
axes[0].set_ylabel('RRI (ms)')
axes[0].set_title('RRI, HR, RMSSD, and Speed Timeseries')
axes[0].grid(True, alpha=0.3)

axes[1].plot(hr_timestamps, hr_values, 'green', linewidth=1)
axes[1].set_ylabel('HR (bpm)')
axes[1].grid(True, alpha=0.3)

axes[2].plot(rmssd_timestamps, rmssd_values, 'orange', linewidth=1)
axes[2].set_ylabel('RMSSD (ms)')
axes[2].grid(True, alpha=0.3)

axes[3].plot(speed_timestamps, speed_values, 'purple', linewidth=1)
axes[3].set_ylabel('Speed')
axes[3].set_xlabel('Time')
axes[3].grid(True, alpha=0.3)

axes[3].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
fig2.autofmt_xdate()
plt.tight_layout()
plt.savefig('rri_hr_rmssd_speed_plot.png', dpi=150)
plt.show()
