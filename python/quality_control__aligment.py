import sys
import dask.dataframe as dd
import matplotlib.pyplot as plt

def main(file_name):
    coverage_file = f"data/qualityControl/alignment/{file_name}.txt"
    # Specify the dtype for the "coverage" column
    dtype = {'coverage': 'float64'}
    coverage_data = dd.read_csv(coverage_file, sep="\t", names=["chromosome", "position", "coverage"], dtype=dtype)

    summary_per_chromosome = coverage_data.groupby("chromosome").agg(
        average_coverage=("coverage", "mean"),
        median_coverage=("coverage", "median"),
        min_coverage=("coverage", "min"),
        max_coverage=("coverage", "max"),
        shuffle="tasks"
    ).compute()

    thresholds = [0, 5, 10, 15]
    
    for threshold in thresholds:
        coverage_above_threshold = (coverage_data["coverage"] > threshold).astype(int)
        over_threshold_percentage = coverage_above_threshold.groupby(coverage_data["chromosome"]).mean() * 100
        summary_per_chromosome[f"over_{threshold}_percent"] = over_threshold_percentage.compute()

    summary_per_chromosome.reset_index(inplace=True)

    summary_csv_file = f"data/qualityControl/alignment/{file_name}_summary_per_chromosome.csv"
    summary_per_chromosome.to_csv(summary_csv_file, index=False)

    plt.bar(summary_per_chromosome["chromosome"], summary_per_chromosome["average_coverage"])
    plt.xlabel("Chromosome")
    plt.ylabel("Average Coverage")
    plt.title("Average Coverage per Chromosome")
    plt.xticks(rotation=45)
    plt.tight_layout()

    plot_image_file = f"data/qualityControl/alignment/{file_name}_coverage_plot.png"
    plt.savefig(plot_image_file)

    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <file_name>")
    else:
        file_name = sys.argv[1]
        main(file_name)
