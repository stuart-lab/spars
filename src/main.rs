use clap::{Parser, Subcommand, ArgGroup};
use std::error::Error;

mod stats;
mod subset;
mod io_utils;

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
    },
    /// Subset the MTX file based on specified rows and columns
    #[command(group(
        ArgGroup::new("indices")
            .required(true)
            .args(&["rows", "cols"]),
    ))]
    Subset {
        /// Input MTX file or directory containing matrix.mtx.gz, features.tsv.gz, and barcodes.tsv
        #[arg(short = 'i', long = "input")]
        input: Option<String>,

        /// Output MTX file or directory name (if directory is used, matrix.mtx.gz, features.tsv.gz, and barcodes.tsv will be subsetted)
        #[arg(short, long)]
        output: String,

        /// File containing row indices to retain (one per line)
        #[arg(long)]
        rows: Option<String>,

        /// File containing column indices to retain (one per line)
        #[arg(long)]
        cols: Option<String>,

        /// Do not reindex the output matrix (keep original indices)
        #[arg(long)]
        no_reindex: bool,
    },
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Stats { input, output_prefix } => {
            stats::compute_stats(&input, &output_prefix)?;
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
    }

    Ok(())
}