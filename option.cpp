#include "ezOptionParser.h"
#include "snpspec.h"

using namespace ez;

void Usage(ezOptionParser & opt)
{
    std::string usage;
    opt.getUsage(usage);
    std::cout << usage;
};

int main(int argc, const char * argv[])
{
    ezOptionParser opt;

    opt.overview = "SNPspec\n=======";
    opt.syntax = "snpspec [OPTIONS]";
    opt.example =
        "snpspec --snps file.txt \n"
        "        --expression file.gct.gz\n"
        "        --null-snps file.txt\n"
        "        --snp-intervals file.bed.gz\n"
        "        --gene-intervals file.bed.gz\n"
        "        --condition file.txt\n"
        "        --out folder/\n"
        "        --slop 250e3\n"
        "        --threads 2\n"
        "        --null-snpsets 100\n"
        "        --min-observations 25\n"
        "        --max-iterations 1e6\n\n";
    opt.footer =
        "SNPspec 0.1  Copyright (C) 2013 Kamil Slowikowski\n"
        "This program is free and without warranty.\n";

    opt.add(
        "", // Default.
        0, // Required?
        0, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Display usage instructions.", // Help description.
        "-h",    // Flag token.
        "--help" // Flag token.
    );

    opt.add(
        "", // Default.
        0, // Required?
        0, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Display version and exit.", // Help description.
        "-v",       // Flag token.
        "--version" // Flag token.
    );

    opt.add(
        "", // Default.
        1, // Required?
        -1, // Number of args expected.
        ',', // Delimiter if expecting multiple args.
        "List of SNPs to test.", // Help description.
        "--snps" // Flag token.
    );

    opt.add(
        "", // Default.
        1, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Gene expression file in GCT format. The Name column "
        "must contain Entrez IDs for genes. These IDs are "
        "looked up in the --genes file.",
        "--expression" // Flag token.
    );

    opt.add(
        "", // Default.
        1, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "BED file with all known gene intervals and two name columns:\n"
        "   - Entrez ID\n"
        "   - HGNC Symbol",
        "--gene-intervals" // Flag token.
    );

    opt.add(
        "", // Default.
        1, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "BED file with all known SNP intervals.",
        "--snp-intervals" // Flag token.
    );

    opt.add(
        "", // Default.
        1, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Names of SNPs to sample when generating random SNP sets "
        "matched on the number of implicated genes. These must be "
        "a subset of the --snp-intervals SNPs.",
        "--null-snps" // Flag token.
    );

    opt.add(
        "", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "List of columns in --expression to condition on before calculating "
        "p-values. Each column in --expression is projected onto each column "
        "listed in this file and its projection is subtracted.",
        "--condition" // Flag token.
    );

    opt.add(
        "", // Default.
        1, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Create output files in this directory.", // Help description.
        "--out" // Flag token.
    );

    opt.add(
        "250000", // Default.
        1, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "If a SNP overlaps no genes, extend its interval this many "
        "nucleotides further and try again. [default: 250000]",
        "--slop" // Flag token.
    );

    opt.add(
        "1", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Number of threads to use. [default: 1]",
        "--threads" // Flag token.
    );

    opt.add(
        "10", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Test this many null matched SNP sets, so you can see how the "
        "distributions of null results for each expression column. "
        "[default: 10]",
        "--null-snpsets" // Flag token.
    );

    opt.add(
        "25", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Stop testing a column in --expression after observing this "
        "many null SNP sets with >= specificity scores. [default: 25]",
        "--min-observations" // Flag token.
    );

    opt.add(
        "1000", // Default.
        0, // Required?
        1, // Number of args expected.
        0, // Delimiter if expecting multiple args.
        "Maximum number of null SNP sets generated. [default: 1000]",
        "--max-iterations" // Flag token.
    );

    // Read the options.
    opt.parse(argc, argv);

    if (opt.isSet("-h")) {
        Usage(opt);
        return 1;
    }

    std::vector<std::string> badOptions;
    int i;
    if (!opt.gotRequired(badOptions)) {
        for (i = 0; i < badOptions.size(); ++i) {
            std::cerr
                    << "ERROR: Missing required option "
                    << badOptions[i] << ".\n\n";
        }
        Usage(opt);
        return 1;
    }

    if (!opt.gotExpected(badOptions)) {
        for (i = 0; i < badOptions.size(); ++i) {
            std::cerr
                    << "ERROR: Got unexpected number of arguments for option "
                    << badOptions[i] << ".\n\n";
        }
        Usage(opt);
        return 1;
    }

    std::vector<std::string>
    user_snpset_files;

    std::string
    expression_file,
    gene_intervals_file,
    snp_intervals_file,
    null_snps_file,
    condition_file,
    out_folder;

    opt.get("--snps")->getStrings(user_snpset_files);
    opt.get("--expression")->getString(expression_file);
    opt.get("--gene-intervals")->getString(gene_intervals_file);
    opt.get("--snp-intervals")->getString(snp_intervals_file);
    opt.get("--null-snps")->getString(null_snps_file);
    opt.get("--condition")->getString(condition_file);
    opt.get("--out")->getString(out_folder);

    // Ensure the files exist.
    for (auto f : user_snpset_files) {
        file_exists(f);
    }
    file_exists(expression_file);
    file_exists(gene_intervals_file);
    file_exists(snp_intervals_file);
    file_exists(null_snps_file);
    // Optional.
    if (condition_file.length() > 0) {
        file_exists(condition_file);
    }

    // Create the output directory.
    mkpath(out_folder);

    int
    threads,
    slop;
    
    long
    null_snpset_replicates,
    min_observations,
    max_iterations;

    opt.get("--threads")->getInt(threads);
    opt.get("--null-snpsets")->getLong(null_snpset_replicates);
    opt.get("--min-observations")->getLong(min_observations);

    double
    slop_d,
    max_iterations_d;

    // Read double so we can pass things like "1e6" and "250e3".
    opt.get("--slop")->getDouble(slop_d);
    opt.get("--max-iterations")->getDouble(max_iterations_d);

    slop = slop_d;
    max_iterations = max_iterations_d;

    if (max_iterations <= 0) {
        std::cerr << "ERROR: Invalid option: --max-iterations " 
                  << max_iterations << std::endl
                  << "This option cannot exceed 1e18.\n";
        exit(EXIT_FAILURE);
    }

    if (min_observations >= max_iterations || min_observations <= 0) {
        std::cerr << "ERROR: Invalid option: --min-observations " 
                  << min_observations << std::endl;
        exit(EXIT_FAILURE);
    }

    // Export all of the options used.
    //opt.exportFile(out_folder + "/args.txt", true);

    // Run the analysis.
    snpspec(
        user_snpset_files,
        expression_file,
        gene_intervals_file,
        snp_intervals_file,
        null_snps_file,
        condition_file,
        out_folder,
        slop,
        threads,
        null_snpset_replicates,
        min_observations,
        max_iterations
    );

    return 0;
}
