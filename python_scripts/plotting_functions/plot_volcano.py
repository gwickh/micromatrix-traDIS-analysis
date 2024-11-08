import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize

#volcano plot
def plot_volcano(
    data, 
    col_name, 
    x_var, 
    y_var, 
    hue_var, 
    style_var, 
    top_n=10,
    label_var='gene_name',
    title="",
    palette={'Downregulated': 'red', 'Upregulated': 'blue', 'Neutral': 'gray'},
    markers={"5'": '^', "3'": 'v', 'Coding': 'X'},
    hue_order=['Downregulated', 'Neutral', 'Upregulated'],  # Order for hue
    style_order=["5'", "3'", 'Coding'],  # Order for style
    s=100, 
    alpha=0.4,
    legend_title='Insertion Position',
    min_dist=0.15,  # Minimum distance between the annotation and the data point
    min_inter_annotation_dist=0.25  # Minimum distance between annotations
):
    # Sort the column values alphabetically for ordering
    col_order = sorted(data[col_name].unique())
    
    # Create FacetGrid with alphabetical ordering
    g = sns.FacetGrid(data, col=col_name, height=6, aspect=0.8, col_order=col_order)
    
    # Map the scatterplot, set zorder for points to be higher than annotations
    g.map(
        sns.scatterplot,
        x_var, 
        y_var,
        legend='full',
        hue=hue_var,
        palette=palette,
        hue_order=hue_order,
        style=style_var,
        markers=markers,
        style_order=style_order,
        s=s,
        alpha=alpha,
        data=data,
        zorder=2  # Data points on top of annotations
    )
    
    # Set custom legend for insertion position
    handles, labels = plt.gca().get_legend_handles_labels()
    # Filter and reorder handles based on style_order for correct legend
    filtered_handles = [handle for handle, label in zip(handles, labels) if label in style_order]
    reordered_handles = [filtered_handles[style_order.index(label)] for label in style_order]

    plt.legend(
        title=legend_title,
        handles=reordered_handles, 
        labels=style_order
    )

    # Set axis labels
    g.set_axis_labels(
        r'log$_2$ Fold Change', 
        r'-log$_{10}$(q-value)'
    )
    
    # Add reference lines and borders
    for ax in g.axes.flat:
        ax.axhline(y=-np.log10(0.05), color='black')
        ax.axvline(x=1, color='black')
        ax.axvline(x=-1, color='black')
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(1)
            spine.set_color('black')
    
    # Set facet titles and main title
    g.set_titles("{col_name}")
    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(title, fontsize=16)
    
    # Label top points within each facet
    for ax in g.axes.flat:
        facet_data = data[data[col_name] == ax.get_title()]
        top_points = facet_data.nlargest(top_n, y_var)

        # Initial positions for annotations
        labels = []
        initial_positions = []

        # Collect data for each top point
        for _, row in top_points.iterrows():
            x, y = row[x_var], row[y_var]
            label = row[label_var]
            labels.append((x, y, label))
            initial_positions.append((x + 0.1, y + 0.1))  # Initial position offset for labels

        # Step 1: Apply the minimum distance constraint from data points
        adjusted_positions = []
        for i, (x, y, label) in enumerate(labels):
            offset_x, offset_y = initial_positions[i]
            # Ensure each annotation is at least `min_dist` away from the point
            while np.sqrt((offset_x - x)**2 + (offset_y - y)**2) < min_dist:
                offset_x += np.random.choice([-1, 1]) * 0.1  # Random step in x direction
                offset_y += np.random.choice([-1, 1]) * 0.1  # Random step in y direction
            adjusted_positions.append((offset_x, offset_y))

        # Step 2: Avoid annotation overlap
        def distance(p1, p2):
            return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

        # Check and adjust positions to maintain minimum distance between annotations
        for i, (x1, y1) in enumerate(adjusted_positions):
            for j, (x2, y2) in enumerate(adjusted_positions):
                if i != j:
                    while distance((x1, y1), (x2, y2)) < min_inter_annotation_dist:
                        # Apply random displacement to avoid overlap
                        x1 += np.random.choice([-1, 1]) * 0.1
                        y1 += np.random.choice([-1, 1]) * 0.1
                        adjusted_positions[i] = (x1, y1)

        # Step 3: Add annotations at adjusted positions (set zorder for annotations to be lower than points)
        for i, (x, y, label) in enumerate(labels):
            optimal_x, optimal_y = adjusted_positions[i]
            ax.text(
                optimal_x, optimal_y, label,
                color='black', fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white", alpha=0.7),
                zorder=1  # Annotations below the data points
            )
    plt.show()