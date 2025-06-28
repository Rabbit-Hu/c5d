import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

key_ratios = [0.90, 0.95, 0.99]
anno_pad = 5

names_mapping = {
    "c5d_linear": "C\\textsuperscript{5}D-Linear",
    "c5d_quad": "C\\textsuperscript{5}D-Quad",
    "c5d_quad_pw": "C\\textsuperscript{5}D-Quad-Pw",
}


def plot_curve(
    ax,
    key_ratio,
    csv_file,
    n_raw,
    has_y_label=False,
    has_legend=False,
    has_title=False,
    has_xlabel=False,
):

    df = pd.read_csv(csv_file)
    func_names = df.columns[1:]
    n_iters = df.shape[0]

    # Qualitative color scheme
    # colors = plt.get_cmap("tab10")
    colors = ["#F87167", "#B260FF", "#60B6F0"]

    ax.set_yscale("log")

    for func_i, func_name in enumerate(func_names):
        ax.plot(df[func_name], label=names_mapping[func_name], color=colors[func_i])
        # mean_iter = np.sum(df[func_name]) / df[func_name][0]
        # ax.axvline(mean_iter, color="black", linestyle="--", alpha=0.5)
        # ax.grid(color="gray", linestyle="--", linewidth=0.5, axis="y")
        ax.grid(color="white", alpha=1, axis="y", linestyle="-", linewidth=1.5)
        ax.set_facecolor((0.94, 0.94, 0.94))

        # despine
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        # set y ticks to 10^x
        plt.yticks([10**i for i in range(0, 4)])
        ax.tick_params(axis="y", which="both", length=0)

    # if has_title:
    #     ax.set_title("$N_\\mathrm{raw}$=" + str(n_raw))

    # if has_xlabel:
    ax.set_xlabel("Iterations")
    if has_y_label:
        ax.set_ylabel("Remaining")
        ax.annotate(
            # "$N_\\mathrm{raw}$=" + str(n_raw),
            "$P_k$=" + str(key_ratio),
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad - anno_pad, 0),
            xycoords=ax.yaxis.label,
            textcoords="offset points",
            size="large",
            ha="right",
            va="center",
            # rotate by 90 degrees
            rotation=90,
        )
    if has_legend:
        leg = ax.legend()
        leg.get_frame().set_linewidth(0.0)


fig = plt.figure(figsize=(9, 2), dpi=300)
plt.rcParams.update(
    {
        "text.usetex": True,
        # Siggraph font libertine
        "font.family": "libertine",
    }
)

exp_root = "outputs/base"
for id, key_ratio in enumerate(key_ratios):
    ax = fig.add_subplot(1, 3, id + 1)
    csv_file = os.path.join(exp_root, "curve_" + str(int(key_ratio * 100)) + ".csv")
    # print(csv_file)
    plot_curve(
        ax,
        key_ratio,
        csv_file,
        10,
        has_y_label=True,
        has_legend=False,
        has_title=id == 0,
        has_xlabel=id==2,
    )

plt.tight_layout()
plt.savefig("curves.pdf")