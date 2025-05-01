use clap::{Parser, Subcommand};
use std::error::Error;
use std::path::PathBuf;
mod stats;
mod subset;
mod io_utils;
mod man;

#[derive(Parser)]
#[command(
    name = env!("CARGO_PKG_NAME"),
    author = env!("CARGO_PKG_AUTHORS"),
    about = env!("CARGO_PKG_DESCRIPTION"),
    version = env!("CARGO_PKG_VERSION")
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compute row and column statistics of the MTX file
    Stats {
        /// Input MTX file
        #[arg(short, long)]
        input: String,

        /// Output prefix for the statistics files
        #[arg(short, long, default_value = "stats")]
        output_prefix: String,

        /// Column to sort the row and column statistics by
        /// Can be one of: NonZeroCount, Sum, Mean, Variance, StdDev, Min, Max, PearsonResidualVar
        #[arg(short, long, default_value = "PearsonResidualVar")]
        sort_by: Option<String>,

        /// Dispersion parameter for Pearson residual variance calculation
        #[arg(short, long, default_value = "100")]
        theta: Option<f64>,
    },
    /// Subset the MTX file based on specified rows and columns
    Subset {
        /// Input MTX file or directory containing matrix.mtx.gz, features.tsv.gz, and barcodes.tsv
        #[arg(short = 'i', long = "input")]
        input: Option<String>,

        /// Output MTX file or directory name (if directory is used, matrix.mtx.gz, features.tsv.gz, and barcodes.tsv will be subsetted)
        #[arg(short, long)]
        output: String,

        /// File containing row indices to retain (one per line)
        #[arg(long, required_unless_present = "cols")]
        rows: Option<String>,

        /// File containing column indices to retain (one per line)
        #[arg(long, required_unless_present = "rows")]
        cols: Option<String>,

        /// Do not reindex the output matrix (keep original indices)
        #[arg(long)]
        no_reindex: bool,
    },

    #[command(
        name = "generate-manpages",
        about = "Generate man pages",
        hide = true  // Hide from normal help output
    )]
    GenerateManPages {
        #[arg(
            short, 
            long,
            help = "Output directory for man pages"
        )]
        outdir: PathBuf,
    },
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Stats { input, output_prefix, sort_by, theta } => {
            stats::compute_stats(&input, &output_prefix, sort_by, theta)?;
        }
        Commands::Subset {
            input,
            output,
            rows,
            cols,
            no_reindex,
        } => {

            subset::subset_matrix(
                input,
                &output,
                rows,
                cols,
                no_reindex,
            )?;
        }
        Commands::GenerateManPages { outdir } => {
            man::generate_manpages(outdir)?;
        },
    }

    Ok(())
}