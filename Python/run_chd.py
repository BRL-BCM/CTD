#!/usr/bin/env python3
"""
CHD Gene Subgroup CTD Analysis Runner
======================================

This module automates the execution of CTD (Comparative Toxicogenomics Database)
analysis on different combinations of CHD gene subgroups, evaluating how gene set
unions improve the identification of connected nodes in disease networks.

Author: Vladimir Kovacevic
Date: 2025
Version: 1.0

Purpose:
--------
Execute CTD_union.py on multiple combinations of CHD gene subgroups to:
1. Identify disease module connectivity patterns
2. Evaluate synergistic effects of combining gene sets
3. Rank combinations by their union improvement metric

Methodology:
------------
For each pair of gene subgroups (S1, S2):
- Run CTD analysis on S1 alone → get F1 nodes
- Run CTD analysis on S2 alone → get F2 nodes
- Run CTD analysis on S1 ∪ S2 → get F_union nodes
- Calculate improvement: |F_union| - max(|F1|, |F2|)

The improvement metric quantifies how much additional disease-relevant connectivity
is discovered by combining gene sets versus using either set individually.

Output:
-------
chd_report.csv with columns:
- combination_id: Unique identifier for the pair
- group1: First gene subgroup filename
- group2: Second gene subgroup filename
- genes_s1: Number of genes in S1
- genes_s2: Number of genes in S2
- genes_union: Number of genes in S1 ∪ S2
- f_nodes_s1: Number of F nodes from S1
- f_nodes_s2: Number of F nodes from S2
- f_nodes_union: Number of F nodes from union
- p_value_s1: Statistical significance for S1
- p_value_s2: Statistical significance for S2
- p_value_union: Statistical significance for union
- improvement: |F_union| - max(|F1|, |F2|)
- improvement_pct: (improvement / max(|F1|, |F2|)) * 100
- penalized_score: -log(p_union) * (|F_union| / max(|F1|, |F2|))

Rows are sorted by improvement (descending) to show most synergistic combinations.

Logging:
--------
The script provides comprehensive logging at INFO and DEBUG levels:
- INFO: Combination selection, CTD execution progress, result collection
- DEBUG: Detailed gene counts, p-values, F node lists

Each step includes:
  ✓ Success indicators
  → Processing indicators
  ■ Result indicators
"""

import argparse
import csv
import itertools
import json
import logging
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Set

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class CTDAnalysisRunner:
    """
    Orchestrates CTD analysis on combinations of CHD gene subgroups.

    This class manages the workflow of:
    1. Selecting gene subgroup combinations
    2. Running CTD_union.py via subprocess
    3. Parsing and aggregating results
    4. Computing improvement metrics
    5. Generating sorted reports
    """

    def __init__(self, adj_matrix: str, gene_groups_dir: str, num_combinations: int = 3):
        """
        Initialize the CTD analysis runner.

        Parameters:
        -----------
        adj_matrix : str
            Path to adjacency matrix CSV file (e.g., adj_TOF.csv)
        gene_groups_dir : str
            Directory containing chd_genes_*.csv files
        num_combinations : int
            Number of gene subgroup pairs to analyze (default: 3)
        """
        self.adj_matrix = Path(adj_matrix)
        self.gene_groups_dir = Path(gene_groups_dir)
        self.num_combinations = num_combinations
        self.results = []

        # Validate inputs
        if not self.adj_matrix.exists():
            raise FileNotFoundError(f"Adjacency matrix not found: {adj_matrix}")
        if not self.gene_groups_dir.exists():
            raise FileNotFoundError(f"Gene groups directory not found: {gene_groups_dir}")

        logger.info(f"Initialized CTD Analysis Runner")
        logger.info(f"  Adjacency matrix: {self.adj_matrix}")
        logger.info(f"  Gene groups dir:  {self.gene_groups_dir}")
        logger.info(f"  Combinations:     {self.num_combinations}")

    def get_gene_group_files(self) -> List[Path]:
        """
        Get all CHD gene subgroup CSV files from the specified directory.

        Returns:
        --------
        List[Path]: Sorted list of gene group file paths

        Notes:
        ------
        Files are sorted to ensure reproducible combination selection.
        Only files matching pattern 'chd_genes_*.csv' are included.
        """
        logger.info("Scanning for gene subgroup files...")

        pattern = "chd_genes_*.csv"
        gene_files = sorted(self.gene_groups_dir.glob(pattern))

        if not gene_files:
            raise FileNotFoundError(
                f"No gene group files found matching pattern '{pattern}' "
                f"in {self.gene_groups_dir}"
            )

        logger.info(f"✓ Found {len(gene_files)} gene subgroup files")
        logger.debug(f"  Files: {[f.name for f in gene_files[:5]]}...")

        return gene_files

    def select_combinations(self, gene_files: List[Path]) -> List[Tuple[Path, Path]]:
        """
        Select diverse combinations of gene subgroups for analysis.

        Strategy:
        ---------
        To maximize biological diversity, we select combinations from different
        categories rather than splitting within the same category:

        1. Group files by category prefix (e.g., 'septal_defects', 'outflow_tract')
        2. Prefer inter-category combinations over intra-category
        3. Select evenly distributed pairs across the category space

        Parameters:
        -----------
        gene_files : List[Path]
            List of all available gene group files

        Returns:
        --------
        List[Tuple[Path, Path]]: Selected combinations as file path pairs

        Notes:
        ------
        If num_combinations exceeds available diverse pairs, falls back to
        including some intra-category combinations.
        """
        logger.info(f"Selecting {self.num_combinations} diverse gene subgroup combinations...")

        # Group files by category (prefix before _N.csv or .csv)
        categories = defaultdict(list)
        for f in gene_files:
            # Extract category: chd_genes_{category}_{N}.csv or chd_genes_{category}.csv
            name = f.stem  # Remove .csv
            parts = name.replace('chd_genes_', '').rsplit('_', 1)

            if len(parts) == 2 and parts[1].isdigit():
                category = parts[0]  # e.g., 'septal_defects' from 'septal_defects_1'
            else:
                category = parts[0]  # e.g., 'general' from 'general'

            categories[category].append(f)

        logger.info(f"  Identified {len(categories)} distinct gene categories")
        logger.debug(f"  Categories: {list(categories.keys())}")

        # Select diverse inter-category combinations
        combinations = []
        category_names = list(categories.keys())

        # Strategy 1: Combinations from different categories (most diverse)
        for i, cat1 in enumerate(category_names):
            for cat2 in category_names[i+1:]:
                # Pick one representative from each category
                file1 = categories[cat1][0]
                file2 = categories[cat2][0]
                combinations.append((file1, file2))

                if len(combinations) >= self.num_combinations:
                    break
            if len(combinations) >= self.num_combinations:
                break

        # Strategy 2: If not enough, add combinations within categories (parts of same category)
        if len(combinations) < self.num_combinations:
            logger.debug("  Adding intra-category combinations to reach target number")
            for cat_name, files in categories.items():
                if len(files) >= 2:
                    for pair in itertools.combinations(files[:3], 2):  # Limit to first 3 parts
                        combinations.append(pair)
                        if len(combinations) >= self.num_combinations:
                            break
                if len(combinations) >= self.num_combinations:
                    break

        # Take only requested number
        selected = combinations[:self.num_combinations]

        logger.info(f"✓ Selected {len(selected)} combinations:")
        for i, (f1, f2) in enumerate(selected, 1):
            logger.info(f"  [{i}] {f1.name} + {f2.name}")

        return selected

    def read_genes_from_csv(self, filepath: Path) -> Set[str]:
        """
        Read gene symbols from a single-column CSV file.

        Parameters:
        -----------
        filepath : Path
            Path to CSV file with column 'x' containing gene symbols

        Returns:
        --------
        Set[str]: Set of gene symbols (uppercase for consistency)

        Notes:
        ------
        - Skips header row (assumes first row is 'x')
        - Converts all genes to uppercase for standardization
        - Ignores empty lines and whitespace
        """
        genes = set()
        with open(filepath, 'r') as f:
            # Skip header
            next(f, None)

            for line in f:
                gene = line.strip()
                if gene and gene.lower() != 'x':
                    genes.add(gene.upper())

        return genes

    def run_ctd_union(
        self,
        s_module1: Path,
        s_module2: Path
    ) -> Dict[str, any]:
        """
        Execute CTD_union.py on a pair of gene subgroup files.

        This runs the CTD analysis three times:
        1. On S_module1 alone
        2. On S_module2 alone
        3. On the union of S_module1 and S_module2

        Parameters:
        -----------
        s_module1 : Path
            First gene subgroup CSV file
        s_module2 : Path
            Second gene subgroup CSV file

        Returns:
        --------
        Dict containing:
            - s_module1_name: Name of first file
            - s_module2_name: Name of second file
            - genes_s1: Number of genes in S1
            - genes_s2: Number of genes in S2
            - genes_union: Number of genes in union
            - results_s1: CTD results for S1 (dict with p_value, F_nodes, etc.)
            - results_s2: CTD results for S2
            - results_union: CTD results for union

        Notes:
        ------
        CTD_union.py creates temporary JSON output files which are parsed
        and included in the returned dictionary.
        """
        logger.info(f"→ Running CTD_union.py on combination:")
        logger.info(f"    S1: {s_module1.name}")
        logger.info(f"    S2: {s_module2.name}")

        # Read gene sets to get sizes
        genes_s1 = self.read_genes_from_csv(s_module1)
        genes_s2 = self.read_genes_from_csv(s_module2)
        genes_union = genes_s1.union(genes_s2)

        logger.debug(f"  Gene counts: |S1|={len(genes_s1)}, |S2|={len(genes_s2)}, "
                    f"|S1∪S2|={len(genes_union)}")

        # Run CTD_union.py
        cmd = [
            'python', 'CTD_union.py',
            '--adj_matrix', str(self.adj_matrix),
            '--s_module', str(s_module1), str(s_module2)
        ]

        logger.debug(f"  Executing: {' '.join(cmd)}")

        try:
            # Capture output
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                cwd=Path.cwd()  # Run from repository root
            )

            logger.debug(f"  CTD_union.py completed successfully")

            # Parse output - CTD_union.py prints results summary
            output_lines = result.stdout.split('\n')

            # Look for JSON output files created by CTD_union.py
            # Pattern: {s_module}.json (without .csv extension, replaced with .json)
            json1 = s_module1.with_suffix('.json')
            json2 = s_module2.with_suffix('.json')

            # Union file is temporary, need to find it from output
            # CTD_union.py creates it as tmpXXX_union.csv and tmpXXX_union.json
            union_json = None
            for line in output_lines:
                if 'union' in line.lower() and '.json' in line:
                    # Extract path from output if possible
                    pass

            # Alternative: parse the printed results
            results = self._parse_ctd_output(result.stdout)

            return {
                's_module1_name': s_module1.name,
                's_module2_name': s_module2.name,
                'genes_s1': len(genes_s1),
                'genes_s2': len(genes_s2),
                'genes_union': len(genes_union),
                'results_s1': results.get(s_module1.name, {}),
                'results_s2': results.get(s_module2.name, {}),
                'results_union': results.get('UNION', {})
            }

        except subprocess.CalledProcessError as e:
            logger.error(f"✗ CTD_union.py failed with exit code {e.returncode}")
            logger.error(f"  STDOUT: {e.stdout}")
            logger.error(f"  STDERR: {e.stderr}")
            raise

    def _parse_ctd_output(self, stdout: str) -> Dict[str, Dict]:
        """
        Parse the summary output from CTD_union.py.

        Expected format:
        === CTD Results Summary ===

        Subset: {filename}
          p-value: {value}
          |F|: {count}
          F nodes (first 10): [gene1, gene2, ...]

        === Penalized P-Value Score ===
          Score: {value}

        Parameters:
        -----------
        stdout : str
            Standard output from CTD_union.py

        Returns:
        --------
        Dict mapping subset name to results dictionary with keys:
            - p_value: float or 'NA'
            - f_nodes_count: int
            - f_nodes_sample: list of first 10 F nodes
            - penalized_score: float or None (for UNION only)
        """
        results = {}
        current_subset = None
        penalized_score = None

        for line in stdout.split('\n'):
            line = line.strip()

            if line.startswith('Subset:'):
                # Extract subset name
                current_subset = line.split('Subset:')[1].strip()
                results[current_subset] = {}

            elif current_subset and line.startswith('p-value:'):
                # Extract p-value
                p_val_str = line.split('p-value:')[1].strip()
                try:
                    p_val = float(p_val_str) if p_val_str != 'NA' else None
                except ValueError:
                    p_val = None
                results[current_subset]['p_value'] = p_val

            elif current_subset and line.startswith('|F|:'):
                # Extract F count
                f_count_str = line.split('|F|:')[1].strip()
                try:
                    f_count = int(f_count_str)
                except ValueError:
                    f_count = 0
                results[current_subset]['f_nodes_count'] = f_count

            elif current_subset and line.startswith('F nodes (first 10):'):
                # Extract F node sample
                nodes_str = line.split('F nodes (first 10):')[1].strip()
                # Parse list format: [gene1, gene2, ...]
                nodes_str = nodes_str.strip('[]')
                if nodes_str:
                    nodes = [n.strip().strip("'\"") for n in nodes_str.split(',')]
                    # Remove trailing '...'
                    nodes = [n for n in nodes if n and n != '...']
                else:
                    nodes = []
                results[current_subset]['f_nodes_sample'] = nodes

            elif line.startswith('Score:') and 'Penalized' in stdout[:stdout.index(line) if line in stdout else 0]:
                # Extract penalized score
                score_str = line.split('Score:')[1].strip()
                if score_str and score_str != 'NA' and '(' not in score_str:
                    try:
                        penalized_score = float(score_str)
                    except ValueError:
                        penalized_score = None

        # Add penalized score to UNION results if available
        if 'UNION' in results and penalized_score is not None:
            results['UNION']['penalized_score'] = penalized_score

        logger.debug(f"  Parsed results for {len(results)} subsets")
        if penalized_score is not None:
            logger.debug(f"  Penalized score: {penalized_score:.4f}")

        return results

    def calculate_improvement(self, result: Dict) -> Dict:
        """
        Calculate improvement metrics for a gene subgroup combination.

        Improvement quantifies the added value of combining gene sets:
        - Absolute improvement: |F_union| - max(|F1|, |F2|)
        - Relative improvement: (absolute / max(|F1|, |F2|)) * 100%

        Positive improvement indicates synergy: the union discovers more
        disease-relevant connections than either set alone.

        Parameters:
        -----------
        result : Dict
            Result dictionary from run_ctd_union()

        Returns:
        --------
        Dict with added keys:
            - f_max_individual: max(|F1|, |F2|)
            - improvement: |F_union| - max(|F1|, |F2|)
            - improvement_pct: (improvement / max(|F1|, |F2|)) * 100
        """
        f1 = result['results_s1'].get('f_nodes_count', 0)
        f2 = result['results_s2'].get('f_nodes_count', 0)
        f_union = result['results_union'].get('f_nodes_count', 0)

        f_max_individual = max(f1, f2)
        improvement = f_union - f_max_individual

        if f_max_individual > 0:
            improvement_pct = (improvement / f_max_individual) * 100
        else:
            improvement_pct = 0.0

        logger.debug(f"  Improvement calculation:")
        logger.debug(f"    |F1| = {f1}, |F2| = {f2}, |F_union| = {f_union}")
        logger.debug(f"    max(|F1|, |F2|) = {f_max_individual}")
        logger.debug(f"    Improvement = {improvement} ({improvement_pct:.1f}%)")

        result['f_max_individual'] = f_max_individual
        result['improvement'] = improvement
        result['improvement_pct'] = improvement_pct

        return result

    def run_analysis(self):
        """
        Execute the full CTD analysis workflow.

        Workflow:
        ---------
        1. Get all available gene subgroup files
        2. Select diverse combinations for analysis
        3. For each combination:
           a. Run CTD_union.py
           b. Parse results
           c. Calculate improvement metrics
        4. Sort results by improvement (descending)
        5. Store in self.results for report generation

        Notes:
        ------
        Progress is logged at each step with clear indicators.
        Errors in individual combinations are logged but don't stop the pipeline.
        """
        logger.info("=" * 70)
        logger.info("Starting CTD Analysis on Gene Subgroup Combinations")
        logger.info("=" * 70)

        # Step 1: Get gene group files
        logger.info("\nStep 1: Loading gene subgroup files...")
        gene_files = self.get_gene_group_files()

        # Step 2: Select combinations
        logger.info("\nStep 2: Selecting gene subgroup combinations...")
        combinations = self.select_combinations(gene_files)

        # Step 3: Run CTD analysis on each combination
        logger.info(f"\nStep 3: Running CTD analysis on {len(combinations)} combinations...")

        for i, (file1, file2) in enumerate(combinations, 1):
            logger.info(f"\n[Combination {i}/{len(combinations)}]")

            try:
                # Run CTD_union.py
                result = self.run_ctd_union(file1, file2)

                # Calculate improvement
                result = self.calculate_improvement(result)

                # Add combination ID
                result['combination_id'] = i

                # Store result
                self.results.append(result)

                logger.info(f"■ Results for combination {i}:")
                logger.info(f"    |F_union| = {result['results_union'].get('f_nodes_count', 0)}")
                logger.info(f"    Improvement = {result['improvement']} "
                          f"({result['improvement_pct']:.1f}%)")

            except Exception as e:
                logger.error(f"✗ Failed to process combination {i}: {str(e)}")
                logger.error(f"  Skipping this combination and continuing...")
                continue

        # Step 4: Sort by improvement
        logger.info("\nStep 4: Sorting results by improvement...")
        self.results.sort(key=lambda x: x['improvement'], reverse=True)

        logger.info(f"✓ Analysis complete. {len(self.results)} combinations processed.")
        logger.info(f"\nTop 3 by improvement:")
        for i, r in enumerate(self.results[:3], 1):
            logger.info(f"  {i}. {r['s_module1_name']} + {r['s_module2_name']}: "
                       f"improvement = {r['improvement']} ({r['improvement_pct']:.1f}%)")

    def generate_report(self, output_file: str = 'chd_report.csv'):
        """
        Generate CSV report of CTD analysis results.

        Report format:
        --------------
        Rows are sorted by improvement (descending) to highlight the most
        synergistic gene subgroup combinations.

        Columns include:
        - Combination identifiers and metadata
        - Gene set sizes
        - F node counts for individual sets and union
        - P-values for statistical significance
        - Improvement metrics (absolute and percentage)

        Parameters:
        -----------
        output_file : str
            Path to output CSV file (default: 'chd_report.csv')

        Notes:
        ------
        P-values are formatted to 6 decimal places or 'NA' if not available.
        All numeric values are properly formatted for downstream analysis.
        """
        logger.info(f"\nStep 5: Generating report to {output_file}...")

        if not self.results:
            logger.warning("⚠ No results to write. Report will be empty.")
            return

        # Define CSV columns
        fieldnames = [
            'combination_id',
            'group1',
            'group2',
            'genes_s1',
            'genes_s2',
            'genes_union',
            'f_nodes_s1',
            'f_nodes_s2',
            'f_nodes_union',
            'p_value_s1',
            'p_value_s2',
            'p_value_union',
            'improvement',
            'improvement_pct',
            'penalized_score'
        ]

        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for result in self.results:
                # Format p-values
                p1 = result['results_s1'].get('p_value')
                p2 = result['results_s2'].get('p_value')
                p_union = result['results_union'].get('p_value')

                def format_pvalue(p):
                    if p is None:
                        return 'NA'
                    try:
                        return f"{float(p):.6f}"
                    except (ValueError, TypeError):
                        return 'NA'

                # Get penalized score from union results
                penalized_score = result['results_union'].get('penalized_score')
                if penalized_score is not None:
                    penalized_score_str = f"{penalized_score:.6f}"
                else:
                    penalized_score_str = 'NA'

                row = {
                    'combination_id': result['combination_id'],
                    'group1': result['s_module1_name'],
                    'group2': result['s_module2_name'],
                    'genes_s1': result['genes_s1'],
                    'genes_s2': result['genes_s2'],
                    'genes_union': result['genes_union'],
                    'f_nodes_s1': result['results_s1'].get('f_nodes_count', 0),
                    'f_nodes_s2': result['results_s2'].get('f_nodes_count', 0),
                    'f_nodes_union': result['results_union'].get('f_nodes_count', 0),
                    'p_value_s1': format_pvalue(p1),
                    'p_value_s2': format_pvalue(p2),
                    'p_value_union': format_pvalue(p_union),
                    'improvement': result['improvement'],
                    'improvement_pct': f"{result['improvement_pct']:.2f}",
                    'penalized_score': penalized_score_str
                }

                writer.writerow(row)

        logger.info(f"✓ Report written successfully")
        logger.info(f"  Rows: {len(self.results)}")
        logger.info(f"  Columns: {len(fieldnames)}")
        logger.info(f"  File: {output_file}")


def main():
    """
    Main execution function for CHD gene subgroup CTD analysis.

    Pipeline:
    ---------
    1. Parse command-line arguments
    2. Initialize CTD analysis runner
    3. Run CTD analysis on selected combinations
    4. Generate sorted report CSV
    5. Display summary statistics

    Exit Codes:
    -----------
    0: Success
    1: Error (file not found, CTD execution failure, etc.)
    """
    parser = argparse.ArgumentParser(
        description='Run CTD analysis on CHD gene subgroup combinations',
        epilog="""
        Example usage:
          %(prog)s --adj_matrix data/scRNA_CHD/adj_TOF.csv --num_combinations 3
          %(prog)s --adj_matrix data/scRNA_CHD/adj_HLHS.csv --num_combinations 5 --output chd_hlhs_report.csv

        Output:
          Creates a CSV report with CTD analysis results sorted by improvement metric.
          The improvement metric quantifies how much the union of two gene sets
          improves disease module connectivity beyond either set alone.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--adj_matrix',
        required=True,
        help='Path to adjacency matrix CSV file (e.g., data/scRNA_CHD/adj_TOF.csv)'
    )
    parser.add_argument(
        '--gene_groups_dir',
        default='data/scRNA_CHD',
        help='Directory containing chd_genes_*.csv files (default: data/scRNA_CHD)'
    )
    parser.add_argument(
        '--num_combinations',
        type=int,
        default=3,
        help='Number of gene subgroup combinations to analyze (default: 3)'
    )
    parser.add_argument(
        '--output',
        default='chd_report.csv',
        help='Output CSV report filename (default: chd_report.csv)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging (DEBUG level)'
    )

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.info("Verbose logging enabled")

    try:
        # Initialize runner
        runner = CTDAnalysisRunner(
            adj_matrix=args.adj_matrix,
            gene_groups_dir=args.gene_groups_dir,
            num_combinations=args.num_combinations
        )

        # Run analysis
        runner.run_analysis()

        # Generate report
        runner.generate_report(output_file=args.output)

        # Print summary
        print(f"\n{'='*78}")
        print(f"{'CTD Analysis Summary':^78}")
        print(f"{'='*78}")
        print(f"  Adjacency matrix:    {args.adj_matrix}")
        print(f"  Combinations tested: {len(runner.results)}")
        print(f"  Report file:         {args.output}")
        print(f"{'='*78}\n")

        if runner.results:
            print("Top 3 combinations by improvement:")
            for i, r in enumerate(runner.results[:3], 1):
                print(f"  {i}. {r['s_module1_name']} + {r['s_module2_name']}")
                print(f"     Improvement: {r['improvement']} nodes ({r['improvement_pct']:.1f}%)")
            print()

        logger.info("=" * 70)
        logger.info("✓ Pipeline completed successfully")
        logger.info("=" * 70)

    except Exception as e:
        logger.error(f"✗ Pipeline failed: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()


"""
IMPLEMENTATION NOTES
====================

CTD Analysis Workflow:
----------------------
The script orchestrates running CTD_union.py, which itself runs CTD.py three times:
1. CTD.py on S1 → identifies F1 (disease-connected nodes for gene set 1)
2. CTD.py on S2 → identifies F2 (disease-connected nodes for gene set 2)
3. CTD.py on S1∪S2 → identifies F_union (disease-connected nodes for combined set)

Improvement Metric:
------------------
improvement = |F_union| - max(|F1|, |F2|)

Interpretation:
- improvement > 0: Synergy - combining sets discovers additional connectivity
- improvement = 0: No synergy - union performs same as best individual set
- improvement < 0: Antagonism - combining sets reduces connectivity (rare)

Biological Significance:
-----------------------
High improvement suggests that the two gene sets target complementary aspects
of the disease network. Genes from S1 and S2 may participate in different
but interconnected pathways that together reveal more of the disease module.

Statistical Considerations:
--------------------------
P-values indicate the statistical significance of the discovered F nodes
representing a true disease module rather than random connectivity. Lower
p-values (typically < 0.05) indicate stronger evidence for disease relevance.
"""
