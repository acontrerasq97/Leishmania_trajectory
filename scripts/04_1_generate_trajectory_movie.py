import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Patch

# 1. Master dictionary mapping base terminal stages to high-contrast colors
terminal_stage_colors = {
    "haptomonads": "#6a3d9a",                # Purple
    "early_metacyclics": "#a6cee3",          # Light Blue
    "late_metacyclics": "#b15928"            # Brown
}
terminal_stages = list(terminal_stage_colors.keys())

# Timeline milestones including Day 5 tracking integration
unique_days_list = [0, 2, 5, 9]

print("-> Initializing High-Definition framework with black-masked background...")
fig, ax = plt.subplots(figsize=(10, 10), facecolor="white")
ax.axis("off")

# Render background layout (full global cellular landscape as a faint context reference)
ax.scatter(
    coord_df["x"],
    coord_df["y"],
    c="#ebdebd", 
    s=3,         
    marker="o",  
    edgecolors="none",
    alpha=0.2,   # Very faint beige backdrop
    zorder=1
)

# Establish the discrete legend panel showcasing both terminal and intermediate categories
legend_elements = [
    Patch(facecolor=color, label=stage.replace("_", " ").title())
    for stage, color in terminal_stage_colors.items()
]
# Append a black patch entry to represent all other developing populations
legend_elements.append(Patch(facecolor="#1a1a1a", label="Intermediate Stages"))

ax.legend(
    handles=legend_elements, 
    loc="lower left", 
    frameon=True, 
    facecolor="white", 
    edgecolor="#e0e0e0",
    title="Leishmania Population Profile",
    fontsize="medium",
    title_fontsize="large"
)

# Initialize foreground container for tracking cell coordinates
foreground_scatter = ax.scatter([], [], s=14, marker="o", edgecolors="none", zorder=3)
title_text = fig.suptitle("Leishmania Terminal Differentiation\nCompiling Video Layers...", fontsize=14, y=0.95)

print("-> Running adaptive chronological mapping loop with black non-terminal styling...")

def animate(i):
    current_day = unique_days_list[i]
    day_suffix = f"_{current_day}"
    
    # Flexible day matching handling float/integer metadata columns smoothly
    try:
        day_mask = np.isclose(trajectory_lmj.obs["day"].astype(float), float(current_day))
    except ValueError:
        day_mask = trajectory_lmj.obs["day"].astype(str) == str(current_day)
        
    _trajectory_subset = trajectory_lmj[day_mask]
    
    if _trajectory_subset.shape[0] == 0:
        title_text.set_text(f"Experimental Timeline: Day {current_day}\nNo Cells Detected")
        foreground_scatter.set_offsets(np.empty((0, 2)))
        return foreground_scatter, title_text
        
    _day_coords = coord_df.loc[_trajectory_subset.obs_names]
    cell_colors = []
    
    # Process categorical cell status mapping across terminal columns exclusively
    for cell_id in _trajectory_subset.obs_names:
        # CHANGED: All intermediate/early cells are now explicitly colored black
        assigned_color = "#1a1a1a"  
        max_val = -float('inf')
        is_terminal = False
        
        for base_stage, color in terminal_stage_colors.items():
            specific_column = f"{base_stage}{day_suffix}"
            
            if specific_column in _trajectory_subset.obs.columns:
                val = _trajectory_subset.obs.loc[cell_id, specific_column]
                
                if isinstance(val, (bool, np.bool_)):
                    if val:
                        assigned_color = color
                        is_terminal = True
                        break
                elif isinstance(val, str):
                    if val.strip().lower() in [base_stage, "true", "1"]:
                        assigned_color = color
                        is_terminal = True
                        break
                else:
                    if val > max_val:
                        max_val = val
                        if val > 0:
                            assigned_color = color
                            is_terminal = True
                            
            elif specific_column in _trajectory_subset.var_names:
                cell_idx = _trajectory_subset.obs_names.get_loc(cell_id)
                stage_idx = _trajectory_subset.var_names.get_loc(specific_column)
                val = _trajectory_subset.X[cell_idx, stage_idx]
                if val > max_val and val > 0:
                    max_val = val
                    assigned_color = color
                    is_terminal = True
        
        # If the cell is active on this day but doesn't map to a terminal score, force black paint
        if not is_terminal:
            assigned_color = "#1a1a1a"
                    
        cell_colors.append(assigned_color)

    # Dynamic title updates tracking terminal census metrics
    title_text.set_text(
        f"Terminal Development: Day {current_day}\n"
        f"Isolating Infective Populations | {_trajectory_subset.shape[0]:,} Total Tracked Cells"
    )
    
    # Refresh foreground visualization layouts
    foreground_scatter.set_offsets(_day_coords[["x", "y"]].values)
    foreground_scatter.set_facecolors(cell_colors)
    
    return foreground_scatter, title_text

print("-> Compiling video frames via FFmpeg...")

anim = animation.FuncAnimation(
    fig,
    func=animate,
    frames=len(unique_days_list),
    blit=False,
    repeat=False
)

movie_file_mp4 = 'leishmania_trajectory/figures/wot/leishmania_terminal_stages_hd.mp4'

# Slow exposure allocation: 3.0 seconds per frame
FFMpegWriter = animation.writers['ffmpeg']
writer = FFMpegWriter(fps=1/3.0, codec='h264', bitrate=3000)

print(f"-> Encoding and exporting HD MP4 to '{movie_file_mp4}'...")
anim.save(movie_file_mp4, writer=writer, dpi=300)

print(f"🎉 Success! Terminal trajectory video with black contrast saved at: '{movie_file_mp4}'")
