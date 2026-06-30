"""Per-mark peak-set Jaccard similarity QC plot.

Computes the pairwise Jaccard index between every peak .bed file for a mark
(each sample's peaks plus the consensus peaks) and renders it as a triangular
heatmap rotated 45 degrees, with a horizontal bar per file showing how many
peaks it contains. Built with plotly.

Each cell is drawn as a real unit square (a filled path) whose corners are
rotated 45 degrees, so neighbouring cells share corners and tile exactly.

Layout (left to right): rotated triangular matrix, peak-count bars, sample names.
"""
import os
import math
import argparse

import numpy as np
import pybedtools
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.colors import sample_colorscale

SQRT2 = math.sqrt(2.0)

# Blues runs light (low Jaccard) -> dark (high Jaccard), which is what we want.
COLORSCALE = "Blues"


def rotate(x, y):
    """Rotate a point 45 degrees so the self-diagonal becomes vertical, with
    the first file at the top (y increases downward, hence the negated y)."""
    return ((x - y) / SQRT2, -(x + y) / SQRT2)


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
            if np.isnan(jaccard):
                jaccard = 0.0
            matrix[row, col] = jaccard
            matrix[col, row] = jaccard
    return matrix


def count_peaks(bed_paths):
    """Number of peaks (data lines) in each .bed file, ignoring headers."""
    counts = []
    for path in bed_paths:
        with open(path) as handle:
            counts.append(sum(
                1 for line in handle
                if line.strip() and not line.startswith(("#", "track", "browser"))
            ))
    return counts


def build_figure(matrix, labels, peak_counts, mark):
    n = len(labels)

    fig = make_subplots(
        rows=1, cols=2, shared_yaxes=True,
        column_widths=[0.5, 0.5], horizontal_spacing=0.06,
    )

    # one filled, rotated unit square per lower-triangle cell (exact tiling)
    centers_x, centers_y, values, hover = [], [], [], []
    for row in range(n):
        for col in range(row + 1):
            value = float(matrix[row, col])
            corners = [(col, row), (col + 1, row), (col + 1, row + 1), (col, row + 1)]
            rotated = [rotate(x, y) for x, y in corners]
            path = "M {} L {} L {} L {} Z".format(
                *(f"{px:.4f},{py:.4f}" for px, py in rotated)
            )
            fig.add_shape(
                type="path", path=path,
                fillcolor=sample_colorscale(COLORSCALE, [value])[0],
                line=dict(color="white", width=0.5),
                row=1, col=1,
            )
            center_x, center_y = rotate(col + 0.5, row + 0.5)
            centers_x.append(center_x)
            centers_y.append(center_y)
            values.append(value)
            hover.append(f"{labels[row]} vs {labels[col]}: {value:.3f}")

    # transparent center markers carry both the colorbar and the hover labels
    fig.add_trace(
        go.Scatter(
            x=centers_x, y=centers_y, mode="markers",
            marker=dict(
                size=12, opacity=0.0,
                color=values, colorscale=COLORSCALE, cmin=0.0, cmax=1.0,
                colorbar=dict(
                    title=dict(text="Jaccard statistic", side="top"),
                    orientation="h", x=0.0, xanchor="left",
                    y=-0.06, yanchor="top", len=0.4, thickness=14,
                ),
                showscale=True,
            ),
            text=hover, hoverinfo="text", showlegend=False,
        ),
        row=1, col=1,
    )

    # diagonal-cell y positions (one per file) for bar + label alignment
    diagonal_y = [rotate(index + 0.5, index + 0.5)[1] for index in range(n)]
    fig.add_trace(
        go.Bar(
            x=peak_counts, y=diagonal_y, orientation="h",
            marker=dict(color="lightskyblue"),
            text=peak_counts, textposition="inside", insidetextanchor="start",
            hovertemplate="%{x} peaks<extra></extra>", showlegend=False,
        ),
        row=1, col=2,
    )

    # matrix axis: keep the rotated squares square, no text labels here
    fig.update_yaxes(
        scaleanchor="x", scaleratio=1,
        showticklabels=False,
        showgrid=False, zeroline=False, row=1, col=1,
    )
    fig.update_xaxes(visible=False, row=1, col=1)

    # bars in the middle; sample names on the right edge of the bar axis
    fig.update_yaxes(
        tickmode="array", tickvals=diagonal_y, ticktext=labels,
        side="right", showgrid=False, zeroline=False, row=1, col=2,
    )
    fig.update_xaxes(
        title_text="Number of peaks", rangemode="tozero",
        showgrid=False, row=1, col=2,
    )

    fig.update_layout(
        title=f"Peak-set Jaccard similarity: {mark}",
        plot_bgcolor="white",
        width=820, height=max(500, 42 * n),
        bargap=0.45, margin=dict(l=20, r=20, t=60, b=110),
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
    # Kaleido has no DPI arg; default is 96 DPI, so scale to reach 300 DPI.
    figure.write_image(args.output_png, scale=300 / 96)
    print(f"Wrote Jaccard QC plot for {args.mark}: {args.output_png}")


if __name__ == "__main__":
    main()