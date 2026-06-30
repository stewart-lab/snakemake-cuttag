"""Per-mark peak-set Jaccard similarity QC plot.

Computes the pairwise Jaccard index between every peak .bed file for a mark
(each sample's peaks plus the consensus peaks) and renders it as a triangular
heatmap rotated 45 degrees, with a horizontal bar per file showing how many
peaks it contains. Built with plotly.
"""
import os
import argparse

import numpy as np
import pybedtools
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def compute_jaccard_matrix(bed_paths):
    """Symmetric n-by-n matrix of pairwise Jaccard indices (diagonal = 1)."""
    sorted_beds = [pybedtools.BedTool(path).sort() for path in bed_paths]
    n = len(sorted_beds)
    matrix = np.ones((n, n))
    for row in range(n):
        for col in range(row):
            try:
                jaccard = sorted_beds[row].jaccard(sorted_beds[col])["jaccard"]
            except Exception:
                # undefined (e.g. an empty peak set) -> treat as no overlap
                jaccard = 0.0
            matrix[row, col] = jaccard
            matrix[col, row] = jaccard
    return matrix


def count_peaks(bed_paths):
    """Number of peaks (non-empty lines) in each .bed file."""
    counts = []
    for path in bed_paths:
        with open(path) as handle:
            counts.append(sum(1 for line in handle if line.strip()))
    return counts


def build_figure(matrix, labels, peak_counts, mark):
    n = len(labels)

    # Rotate the lower triangle 45 degrees: cell (row, col) with row >= col maps
    # to (x = col - row, y = row + col). The self-comparisons (row == col) land
    # on x = 0 at y = 2*row, forming the vertical right edge of the triangle;
    # the triangle fans out to the left.
    xs, ys, values, hover = [], [], [], []
    for row in range(n):
        for col in range(row + 1):
            xs.append(col - row)
            ys.append(row + col)
            values.append(matrix[row, col])
            hover.append(f"{labels[row]} vs {labels[col]}: {matrix[row, col]:.3f}")

    diagonal_y = [2 * index for index in range(n)]

    # marker size (px) is approximate; scale it down as the matrix grows
    marker_size = max(6, int(480 / (n + 1)))

    fig = make_subplots(
        rows=1, cols=2, column_widths=[0.72, 0.28], horizontal_spacing=0.04
    )

    # rotated-square heatmap cells
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode="markers",
            marker=dict(
                symbol="diamond",
                size=marker_size,
                color=values,
                colorscale="RdBu_r",
                cmin=0.0,
                cmax=1.0,
                line=dict(width=0),
                colorbar=dict(
                    title=dict(text="Jaccard statistic", side="top"),
                    orientation="h",
                    x=0.0,
                    xanchor="left",
                    y=-0.02,
                    yanchor="top",
                    len=0.4,
                    thickness=14,
                ),
                showscale=True,
            ),
            text=hover,
            hoverinfo="text",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # one bar per file: number of peaks
    fig.add_trace(
        go.Bar(
            x=peak_counts,
            y=diagonal_y,
            orientation="h",
            marker=dict(color="lightskyblue"),
            text=peak_counts,
            textposition="outside",
            hovertemplate="%{x} peaks<extra></extra>",
            showlegend=False,
        ),
        row=1,
        col=2,
    )

    # share the vertical extent; first file at the top (reversed range)
    y_range = [2 * (n - 1) + 1, -1]
    fig.update_yaxes(
        range=y_range, showticklabels=False, showgrid=False, zeroline=False,
        row=1, col=1,
    )
    fig.update_xaxes(
        showticklabels=False, showgrid=False, zeroline=False, row=1, col=1,
    )
    fig.update_yaxes(
        range=y_range,
        tickmode="array",
        tickvals=diagonal_y,
        ticktext=labels,
        side="left",
        showgrid=False,
        row=1, col=2,
    )
    fig.update_xaxes(
        title_text="Number of peaks", rangemode="tozero", row=1, col=2,
    )

    fig.update_layout(
        title=f"Peak-set Jaccard similarity: {mark}",
        plot_bgcolor="white",
        width=820,
        height=max(420, 55 * n),
        bargap=0.45,
        margin=dict(l=20, r=20, t=60, b=90),
    )
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--beds", nargs="+", required=True,
        help="Peak .bed files for one mark (samples, then consensus).",
    )
    parser.add_argument(
        "--labels", nargs="+", required=True,
        help="Row labels matching the order of --beds.",
    )
    parser.add_argument("--mark", required=True, help="Mark name (plot title).")
    parser.add_argument("--output-png", required=True, help="Output .png path.")
    parser.add_argument("--output-html", required=True, help="Output .html path.")
    args = parser.parse_args()

    if len(args.beds) != len(args.labels):
        parser.error(
            f"--beds ({len(args.beds)}) and --labels ({len(args.labels)}) "
            "must have the same length"
        )

    for output_path in (args.output_png, args.output_html):
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

    matrix = compute_jaccard_matrix(args.beds)
    peak_counts = count_peaks(args.beds)
    figure = build_figure(matrix, args.labels, peak_counts, args.mark)

    figure.write_html(args.output_html)
    figure.write_image(args.output_png)
    print(f"Wrote Jaccard QC plot for {args.mark}: {args.output_png}")


if __name__ == "__main__":
    main()
