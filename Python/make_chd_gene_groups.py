#!/usr/bin/env python3
"""
CHD Gene Subgroup Generator
============================

This module creates biologically and physiologically meaningful subgroups of genes
associated with Congenital Heart Defects (CHD) from scRNA-seq data.

Author: Vladimir Kovacevic
Date: 2025
Version: 1.0

Purpose:
--------
Organize CHD-associated genes into anatomically and functionally coherent groups
to facilitate targeted analysis of disease mechanisms, developmental pathways,
and therapeutic targets.

Input:
------
CSV file containing CHD gene evidence with columns:
- Gene: Gene symbol
- Species: Organism (filtered for 'Human')
- Associated CHD: Semicolon-separated list of CHD types

Output:
-------
Multiple CSV files (one per group) with:
- Filename pattern: chd_genes_{group_name}.csv
- Single column 'x' containing gene symbols
- Maximum of max_genes_per_group genes per file

Logging:
--------
The script provides comprehensive logging at multiple levels:
- INFO level (default): Pipeline steps, category sizes, file creation progress
- DEBUG level (--verbose): Detailed gene assignments, progress percentages,
  species distributions, gene ranges per file

Each major step is numbered and logged with status indicators:
  ✓ Success
  ✗ Failure/Error
  ⚠ Warning
  ✂ Split operation

References:
-----------
CHD classification based on anatomical and developmental considerations
following the International Paediatric and Congenital Cardiac Code (IPCCC).
"""

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class CHDGeneGrouper:
    """
    Organizes CHD-associated genes into biologically meaningful categories.

    This class implements a classification system based on anatomical location,
    developmental origin, and functional characteristics of CHD subtypes.
    """

    # Comprehensive CHD abbreviations dictionary for reference
    CHD_ABBREVIATIONS = {
        # Septal defects
        'VSD': 'Ventricular Septal Defect',
        'ASD': 'Atrial Septal Defect',
        'AVSD': 'Atrioventricular Septal Defect',

        # Conotruncal defects
        'TOF': 'Tetralogy of Fallot',
        'TGA': 'Transposition of the Great Arteries',
        'CTD': 'Conotruncal Defects',
        'DORV': 'Double Outlet Right Ventricle',
        'PTA': 'Persistent Truncus Arteriosus',

        # Pulmonary defects
        'PA': 'Pulmonary Atresia',
        'PS': 'Pulmonary Stenosis',
        'POF': 'Pulmonary Outflow',

        # Aortic defects
        'BAV': 'Bicuspid Aortic Valve',
        'CoA': 'Coarctation of the Aorta',
        'IAA': 'Interrupted Aortic Arch',

        # Left-sided lesions
        'HLHS': 'Hypoplastic Left Heart Syndrome',
        'LVOTO': 'Left Ventricular Outflow Tract Obstruction',
        'LSLs': 'Left-Sided Lesions',

        # Complex defects
        'SV': 'Single Ventricle',
        'ECD': 'Endocardial Cushion Defect',

        # Vascular defects
        'PDA': 'Patent Ductus Arteriosus',
        'RAA': 'Right Aortic Arch',

        # Valve defects
        'MVS': 'Mitral Valve Stenosis',
        'MR': 'Mitral Regurgitation',
        'TV': 'Tricuspid Valve',
        'TA': 'Tricuspid Atresia',

        # Pulmonary venous
        'APVC': 'Anomalous Pulmonary Venous Connection',
        'TAPVC': 'Total Anomalous Pulmonary Venous Connection',
        'PAPVC': 'Partial Anomalous Pulmonary Venous Connection',

        # Laterality
        'HTX': 'Heterotaxy',

        # Other
        'EA': 'Ebstein Anomaly',
        'ALCAPA': 'Anomalous Left Coronary Artery from Pulmonary Artery',
    }

    @staticmethod
    def define_gene_categories() -> Dict[str, Dict[str, any]]:
        """
        Define biologically and physiologically meaningful gene categories.

        Categories are organized based on:
        1. Anatomical location (chambers, valves, vessels)
        2. Developmental origin (endocardial, myocardial, neural crest)
        3. Functional systems (septation, valve formation, laterality)

        Returns:
        --------
        Dict[str, Dict]: Dictionary of categories with metadata
            - description: Human-readable category description
            - biological_rationale: Scientific explanation
            - patterns: CHD type keywords for assignment
            - genes: Set to store assigned genes (initialized empty)
        """

        categories = {
            'septal_defects': {
                'description': 'Septal Defects',
                'biological_rationale': """
                    Genes associated with cardiac septation defects (VSD, ASD, AVSD).

                    DEVELOPMENTAL CONTEXT:
                    Cardiac septation is a complex process occurring during weeks 4-8 of
                    embryonic development, involving the formation of the atrial septum
                    primum/secundum and ventricular septum. Defects arise from incomplete
                    fusion or abnormal remodeling of endocardial cushions and myocardial
                    tissues.

                    MOLECULAR MECHANISMS:
                    - Endocardial-to-mesenchymal transition (EndMT)
                    - Cell migration and proliferation in cushion tissues
                    - Extracellular matrix remodeling
                    - TGF-β, BMP, and Notch signaling pathways

                    CLINICAL SIGNIFICANCE:
                    Most common CHD subtype (~40% of all CHD cases). Can occur in isolation
                    or as part of syndromes (e.g., Down syndrome, 22q11.2 deletion).
                """,
                'patterns': ['VSD', 'ASD', 'AVSD', 'Septal defects'],
                'genes': set()
            },

            'outflow_tract': {
                'description': 'Outflow Tract and Conotruncal Defects',
                'biological_rationale': """
                    Genes involved in outflow tract (OFT) formation and conotruncal defects
                    (TOF, TGA, DORV, CTD).

                    DEVELOPMENTAL CONTEXT:
                    The cardiac OFT develops from the second heart field and undergoes complex
                    remodeling including septation, rotation, and alignment with ventricles.
                    Neural crest cells migrate into the OFT and contribute to aorticopulmonary
                    septation and semilunar valve formation (weeks 5-8).

                    MOLECULAR MECHANISMS:
                    - Neural crest cell migration and differentiation
                    - Second heart field deployment
                    - Retinoic acid, FGF, and Wnt signaling
                    - OFT rotation and wedging mechanisms

                    CLINICAL SIGNIFICANCE:
                    TOF (most common cyanotic CHD) and TGA (medical emergency requiring
                    immediate intervention). Often associated with 22q11.2 deletion syndrome
                    (DiGeorge syndrome).
                """,
                'patterns': ['TOF', 'TGA', 'CTD', 'DORV', 'Malformation of the outflow tracts',
                           'PTA', 'Truncus'],
                'genes': set()
            },

            'pulmonary_defects': {
                'description': 'Pulmonary Valve and Artery Defects',
                'biological_rationale': """
                    Genes affecting pulmonary valve development and right ventricular outflow
                    (PA, PS, pulmonary valve dysplasia).

                    DEVELOPMENTAL CONTEXT:
                    Pulmonary valve develops from endocardial cushions of the truncus
                    arteriosus through condensation, excavation, and remodeling processes.
                    Pulmonary artery formation involves both cardiac neural crest and second
                    heart field contributions.

                    MOLECULAR MECHANISMS:
                    - Semilunar valve morphogenesis and remodeling
                    - Endothelial-to-mesenchymal transition
                    - VEGF signaling in vascular development
                    - Notch pathway in valve stratification

                    CLINICAL SIGNIFICANCE:
                    PS is common (~8-10% of CHD). PA represents severe end of spectrum with
                    complete obstruction. Often requires intervention (valvuloplasty or
                    surgical repair).
                """,
                'patterns': ['PA', 'PS', 'POF', 'Dysplasia of pulmonary valve',
                           'Pulmonary atresia', 'Pulmonary stenosis'],
                'genes': set()
            },

            'aortic_defects': {
                'description': 'Aortic Valve and Arch Defects',
                'biological_rationale': """
                    Genes involved in aortic valve formation and aortic arch remodeling
                    (BAV, CoA, IAA).

                    DEVELOPMENTAL CONTEXT:
                    Aortic valve formation parallels pulmonary valve development. Aortic
                    arch remodeling involves asymmetric regression of pharyngeal arch arteries
                    (weeks 4-7), with persistence of left 4th arch as definitive aortic arch.

                    MOLECULAR MECHANISMS:
                    - Hemodynamic forces shaping valve morphology
                    - Neural crest contribution to aortic arch arteries
                    - Apoptosis in arch artery remodeling
                    - Notch signaling in bicuspid vs. tricuspid valve fate

                    CLINICAL SIGNIFICANCE:
                    BAV is most common congenital heart malformation (~1-2% population).
                    CoA often associated with Turner syndrome and BAV. Can cause heart
                    failure if severe.
                """,
                'patterns': ['BAV', 'CoA', 'IAA', 'Bicuspid aortic', 'Coarctation',
                           'Aortic arch', 'Interrupted aortic'],
                'genes': set()
            },

            'left_sided_lesions': {
                'description': 'Left-Sided Obstructive Lesions',
                'biological_rationale': """
                    Genes causing left heart hypoplasia and obstruction
                    (HLHS, LVOTO, mitral/aortic stenosis).

                    DEVELOPMENTAL CONTEXT:
                    Left-sided lesions represent a spectrum from mild obstruction to complete
                    HLHS. Development involves complex interplay between hemodynamics and
                    genetic factors affecting left ventricle, mitral valve, and aortic valve
                    development.

                    MOLECULAR MECHANISMS:
                    - Reduced blood flow affecting chamber growth (hemodynamic hypothesis)
                    - Primary genetic defects in cardiomyocyte proliferation
                    - Notch, Hippo, and growth factor signaling
                    - Transcription factors regulating chamber specification

                    CLINICAL SIGNIFICANCE:
                    HLHS is severe, requiring staged surgical palliation (Norwood, Glenn,
                    Fontan). Represents ~2-3% of CHD but ~25% of neonatal cardiac deaths.
                    Strong familial clustering suggests genetic basis.
                """,
                'patterns': ['HLHS', 'LVOTO', 'LSLs', 'Left ventricular hypoplasia',
                           'Hypoplastic left', 'Left-sided lesions'],
                'genes': set()
            },

            'single_ventricle': {
                'description': 'Single Ventricle and Complex Defects',
                'biological_rationale': """
                    Genes associated with single ventricle physiology and complex
                    multicomponent CHD (SV, endocardial cushion defects).

                    DEVELOPMENTAL CONTEXT:
                    Single ventricle encompasses heterogeneous defects where one ventricle
                    is absent/hypoplastic or where both atria connect to single ventricle.
                    Represents failure of normal chamber septation and/or atrioventricular
                    valve development.

                    MOLECULAR MECHANISMS:
                    - Defective trabeculation and compaction
                    - Abnormal endocardial cushion formation
                    - Disrupted hand/heart transcription factor networks
                    - Laterality signaling defects (left-right patterning)

                    CLINICAL SIGNIFICANCE:
                    Requires Fontan palliation (single ventricle supporting systemic
                    circulation). Long-term complications include heart failure, arrhythmias,
                    and protein-losing enteropathy.
                """,
                'patterns': ['SV', 'ECD', 'Univentricular', 'Single ventricle',
                           'Endocardial cushion'],
                'genes': set()
            },

            'vascular_defects': {
                'description': 'Patent Ductus and Vascular Anomalies',
                'biological_rationale': """
                    Genes affecting ductus arteriosus closure and major vessel development
                    (PDA, vena cava anomalies, aortic arch variants).

                    DEVELOPMENTAL CONTEXT:
                    Ductus arteriosus is essential for fetal circulation, normally closing
                    within days after birth via smooth muscle contraction and endothelial
                    remodeling. Vena cava and systemic veins derive from sinus venosus and
                    cardinal veins.

                    MOLECULAR MECHANISMS:
                    - Prostaglandin signaling in ductal patency/closure
                    - Oxygen tension sensing and response
                    - Vascular smooth muscle differentiation
                    - VEGF and angiopoietin in vessel remodeling

                    CLINICAL SIGNIFICANCE:
                    PDA common in premature infants (~8 per 1000 births). Can cause heart
                    failure if large. Vena cava anomalies often associated with heterotaxy
                    syndromes.
                """,
                'patterns': ['PDA', 'Vena cava abnormality', 'RAA', 'Patent ductus',
                           'Right aortic arch', 'Vascular ring'],
                'genes': set()
            },

            'valve_abnormalities': {
                'description': 'Atrioventricular Valve Malformations',
                'biological_rationale': """
                    Genes involved in mitral and tricuspid valve development
                    (mitral stenosis/regurgitation, tricuspid dysplasia, Ebstein anomaly).

                    DEVELOPMENTAL CONTEXT:
                    Atrioventricular (AV) valves form from endocardial cushions through
                    EndMT, proliferation, remodeling, and delamination processes (weeks 5-9).
                    Valve leaflets undergo thinning and stratification to form functional
                    tricuspid (right) and bicuspid/mitral (left) valves.

                    MOLECULAR MECHANISMS:
                    - TGF-β and BMP signaling in EndMT
                    - Extracellular matrix organization (collagens, proteoglycans)
                    - Wnt and Notch in valve remodeling
                    - Has2-mediated hyaluronan production

                    CLINICAL SIGNIFICANCE:
                    Ebstein anomaly (apical displacement of tricuspid valve) can range from
                    mild to severe with right heart failure. Mitral valve prolapse affects
                    ~2% population. Often require surgical repair or replacement.
                """,
                'patterns': ['MVS', 'MR', 'TV', 'Congenital anomaly of tricuspid valve',
                           'Dysplasia of mitral valve', 'Ebstein', 'EA',
                           'Mitral valve', 'Tricuspid valve'],
                'genes': set()
            },

            'pulmonary_venous': {
                'description': 'Pulmonary Venous Connection Defects',
                'biological_rationale': """
                    Genes causing abnormal pulmonary venous drainage
                    (TAPVC, PAPVC, anomalous pulmonary venous connections).

                    DEVELOPMENTAL CONTEXT:
                    Pulmonary veins initially drain to systemic venous system. Normal
                    development requires formation of common pulmonary vein, incorporation
                    into left atrium, and regression of embryonic connections to cardinal/
                    umbilicovitelline systems (weeks 4-7).

                    MOLECULAR MECHANISMS:
                    - Pitx2 and other laterality genes in venous patterning
                    - Vascular endothelial guidance cues
                    - Apoptosis in remodeling embryonic connections
                    - Hemodynamic influences on vessel persistence

                    CLINICAL SIGNIFICANCE:
                    TAPVC requires urgent surgical repair (obstructed forms are emergencies).
                    PAPVC often discovered incidentally, may cause right heart volume
                    overload. Can be associated with sinus venosus ASD.
                """,
                'patterns': ['APVC', 'TAPVC', 'PAPVC', 'Anomalous pulmonary venous',
                           'Pulmonary venous'],
                'genes': set()
            },

            'laterality_defects': {
                'description': 'Heterotaxy and Laterality Defects',
                'biological_rationale': """
                    Genes regulating left-right body axis and cardiac situs
                    (heterotaxy, isomerism, situs inversus with CHD).

                    DEVELOPMENTAL CONTEXT:
                    Left-right asymmetry is established early in development through nodal
                    cilia-driven leftward flow, asymmetric gene expression (Nodal, Lefty,
                    Pitx2), and subsequent organ lateralization. Defects cause randomized
                    or abnormal organ positioning.

                    MOLECULAR MECHANISMS:
                    - Nodal cilia structure and motility (dynein arms)
                    - Nodal/TGF-β signaling cascade
                    - Pitx2 transcription factor in left-sided identity
                    - Retinoic acid gradients

                    CLINICAL SIGNIFICANCE:
                    Heterotaxy syndrome involves complex CHD (AVSD, TAPVC, interrupted IVC)
                    plus situs abnormalities of abdominal organs. Associated with asplenia
                    or polysplenia. High morbidity/mortality due to complex defects.
                """,
                'patterns': ['Isomerism', 'Heterotaxy', 'HTX', 'Situs', 'Laterality'],
                'genes': set()
            },

            'right_sided_lesions': {
                'description': 'Right Ventricular and Tricuspid Defects',
                'biological_rationale': """
                    Genes causing right ventricular hypoplasia and tricuspid valve defects
                    (tricuspid atresia, RV hypoplasia, Uhl anomaly).

                    DEVELOPMENTAL CONTEXT:
                    Right-sided lesions mirror left-sided spectrum but are less common.
                    Right ventricle develops from both first and second heart fields.
                    Severe forms (tricuspid atresia) result in single ventricle physiology.

                    MOLECULAR MECHANISMS:
                    - Chamber-specific transcription factors (Irx4, Hand1/2)
                    - Trabeculation and compaction signaling
                    - Valve morphogenesis pathways
                    - Hemodynamic influences on RV growth

                    CLINICAL SIGNIFICANCE:
                    Tricuspid atresia represents ~1-3% of CHD, requires single ventricle
                    palliation. RV hypoplasia often part of pulmonary atresia with intact
                    ventricular septum (PA-IVS) complex.
                """,
                'patterns': ['Right ventricular hypoplasia', 'TA', 'Tricuspid atresia',
                           'RV hypoplasia', 'Uhl'],
                'genes': set()
            },

            'coronary_anomalies': {
                'description': 'Coronary Artery Anomalies',
                'biological_rationale': """
                    Genes affecting coronary artery development and connection
                    (ALCAPA, coronary fistulas, anomalous coronary origins).

                    DEVELOPMENTAL CONTEXT:
                    Coronary arteries develop through vasculogenesis in epicardium, then
                    connect to aortic sinuses via endothelial invasion (weeks 7-9). Process
                    requires precise guidance of coronary stems to appropriate aortic sinuses.

                    MOLECULAR MECHANISMS:
                    - Epicardial-derived cell differentiation
                    - VEGF signaling in coronary angiogenesis
                    - Endothelial guidance cues for aortic connection
                    - Hemodynamic shear stress in vessel remodeling

                    CLINICAL SIGNIFICANCE:
                    ALCAPA causes myocardial ischemia and heart failure in infancy (requires
                    urgent surgery). Anomalous coronary origins from opposite sinus can cause
                    sudden cardiac death in young athletes.
                """,
                'patterns': ['Coronary', 'ALCAPA', 'Anomalous coronary'],
                'genes': set()
            },
        }

        return categories

    @staticmethod
    def parse_chd_types(chd_string: str) -> List[str]:
        """
        Parse and clean CHD type strings from 'Associated CHD' column.

        Parameters:
        -----------
        chd_string : str
            Semicolon-separated string of CHD types

        Returns:
        --------
        List[str]: List of individual CHD type strings

        Examples:
        ---------
        >>> parse_chd_types("VSD; ASD; TOF")
        ['VSD', 'ASD', 'TOF']
        >>> parse_chd_types("CHD")
        ['CHD']
        """
        if pd.isna(chd_string):
            return []

        # Split by semicolon and remove leading/trailing whitespace
        types = [chd_type.strip() for chd_type in str(chd_string).split(';')]
        return types

    def categorize_genes(
        self,
        df: pd.DataFrame,
        max_groups: int = 20,
        max_genes_per_group: int = 50
    ) -> Tuple[Dict, Set[str], Dict[str, List[str]]]:
        """
        Assign genes to biological/physiological categories based on associated CHD.

        Algorithm:
        ----------
        1. For each gene-CHD association in the dataset:
           a. Parse the CHD type(s) from 'Associated CHD' column
           b. Match CHD types against predefined category patterns
           c. Assign gene to all matching categories (genes can be multi-assigned)

        2. Track gene assignments to enable analysis of:
           - Genes involved in multiple developmental processes
           - Shared molecular pathways across CHD types
           - Pleiotropic effects of genetic variants

        Parameters:
        -----------
        df : pd.DataFrame
            Input dataframe with Gene, Species, and Associated CHD columns
        max_groups : int
            Maximum number of output groups (informational)
        max_genes_per_group : int
            Target maximum genes per group (groups may be split)

        Returns:
        --------
        Tuple containing:
            - categories (Dict): Populated category dictionary with gene sets
            - general_chd (Set[str]): Genes with non-specific "CHD" annotation
            - gene_assignments (Dict[str, List[str]]): Gene -> category mappings

        Notes:
        ------
        Genes can be assigned to multiple categories because:
        - A gene may have multiple CHD associations in the literature
        - Same gene may affect multiple developmental processes
        - Pleiotropic genes participate in shared pathways
        """
        logger.info("Initializing gene categorization...")

        # Get predefined categories
        categories = self.define_gene_categories()

        # General CHD category for non-specific associations
        general_chd = set()

        # Track assignments for each gene (for analysis purposes)
        gene_assignments = defaultdict(list)

        # Counters for logging
        total_associations = 0
        categorized_associations = 0

        # Process each gene-CHD association
        logger.info(f"Processing {len(df)} gene-CHD associations...")

        # Progress tracking
        progress_interval = max(1, len(df) // 10)  # Log every 10% of progress

        for idx, row in df.iterrows():
            gene = row['Gene']
            chd_types = self.parse_chd_types(row['Associated CHD'])

            # Log progress at intervals
            if idx > 0 and idx % progress_interval == 0:
                progress_pct = (idx / len(df)) * 100
                logger.debug(f"Progress: {progress_pct:.0f}% ({idx}/{len(df)} entries processed)")

            # Track if this association matched any specific category
            matched_specific_category = False

            for chd_type in chd_types:
                total_associations += 1
                logger.debug(f"  Gene: {gene}, CHD type: {chd_type}")

                # Check against each predefined category
                for cat_name, cat_info in categories.items():
                    for pattern in cat_info['patterns']:
                        # Pattern matching: substring match allows flexibility
                        # e.g., "VSD" matches "VSD", "Perimembranous VSD", etc.
                        if pattern in chd_type:
                            categories[cat_name]['genes'].add(gene)
                            gene_assignments[gene].append(cat_name)
                            matched_specific_category = True
                            categorized_associations += 1
                            logger.debug(f"    → Matched to category: {cat_name} (pattern: {pattern})")
                            break  # Stop after first pattern match per category

                # Assign to general CHD if only non-specific annotation
                if not matched_specific_category and chd_type.strip() == 'CHD':
                    general_chd.add(gene)
                    logger.debug(f"    → Assigned to general CHD")

        # Log categorization statistics
        logger.info(f"Processed {total_associations} gene-CHD associations")
        logger.info(f"Categorized {categorized_associations} specific associations")
        logger.info(f"Found {len(general_chd)} genes with general CHD annotation only")

        # Log category sizes
        for cat_name, cat_info in categories.items():
            if len(cat_info['genes']) > 0:
                logger.info(f"  {cat_name}: {len(cat_info['genes'])} genes")

        return categories, general_chd, gene_assignments

    def create_intersection_groups(
        self,
        categories: Dict,
        max_genes_per_group: int,
        gene_assignments: Dict[str, List[str]]
    ) -> Dict[str, Dict]:
        """
        Create biologically meaningful groups through category intersections.

        Strategy:
        ---------
        Pure intersection-based approach (NO SPLITTING):
        1. For large categories (>max_genes_per_group): Create intersections with related categories
        2. Keep creating intersections until groups are small enough or no more intersections possible
        3. Never split categories - only merge through intersections
        4. Small categories (≤max_genes_per_group): Keep intact

        Biological Rationale:
        --------------------
        Genes appearing in multiple categories represent:
        - Pleiotropic genes affecting multiple developmental processes
        - Shared signaling pathways between cardiac structures
        - Master regulators with broad developmental roles
        - Genes at pathway intersections

        Intersection Combinations:
        -------------------------
        - Septal defects × Outflow tract (shared endocardial cushion origins)
        - Aortic × Left-sided lesions (left heart developmental continuity)
        - Pulmonary × Vascular defects (pulmonary vasculature development)
        - Outflow tract × Pulmonary defects (right ventricular outflow)
        - Outflow tract × Aortic defects (semilunar valve development)
        - Valve abnormalities × Septal defects (endocardial cushion shared origin)
        - Single ventricle × Laterality defects (complex developmental defects)
        - Vascular × Pulmonary venous (systemic/pulmonary venous development)

        Parameters:
        -----------
        categories : Dict
            Category dictionary with gene sets
        max_genes_per_group : int
            Target maximum genes per output file (used as threshold, not hard limit)
        gene_assignments : Dict[str, List[str]]
            Gene-to-category assignments for finding intersections

        Returns:
        --------
        Dict[str, Dict]: Final group dictionary with intersection groups only (no splits)
        """
        logger.info(f"Creating intersection-based groups (target max {max_genes_per_group} genes)...")
        logger.info(f"Strategy: Pure intersections, NO splitting")
        logger.info(f"Analyzing {len(categories)} categories...")

        final_groups = {}

        # Show category sizes
        for cat_name, cat_info in categories.items():
            if len(cat_info['genes']) > 0:
                logger.info(f"  {cat_name}: {len(cat_info['genes'])} genes")

        # Define biologically meaningful intersection pairs
        # Format: (category1, category2, description, rationale)
        intersection_pairs = [
            ('septal_defects', 'outflow_tract',
             'Septal & Outflow Tract Defects',
             'Genes affecting both septation and conotruncal development (shared endocardial cushion origin)'),

            ('aortic_defects', 'left_sided_lesions',
             'Aortic & Left Heart Defects',
             'Genes affecting left ventricular outflow tract and left heart development'),

            ('pulmonary_defects', 'vascular_defects',
             'Pulmonary & Vascular Defects',
             'Genes involved in pulmonary artery and ductus arteriosus development'),

            ('outflow_tract', 'pulmonary_defects',
             'Conotruncal & Pulmonary Defects',
             'Genes affecting right ventricular outflow tract and pulmonary valve'),

            ('outflow_tract', 'aortic_defects',
             'Conotruncal & Aortic Defects',
             'Genes involved in semilunar valve development (aortic and pulmonary)'),

            ('valve_abnormalities', 'septal_defects',
             'Valve & Septal Defects',
             'Genes affecting endocardial cushion-derived structures (AV valves and septa)'),

            ('single_ventricle', 'laterality_defects',
             'Complex Cardiac & Laterality Defects',
             'Genes causing complex cardiac malformations with laterality defects'),

            ('vascular_defects', 'pulmonary_venous',
             'Systemic & Pulmonary Venous Defects',
             'Genes affecting systemic and pulmonary venous development'),

            ('septal_defects', 'valve_abnormalities',
             'Septal & AV Valve Defects',
             'Genes affecting atrioventricular canal development (AVSD spectrum)'),

            ('left_sided_lesions', 'valve_abnormalities',
             'Left Heart & Valve Defects',
             'Genes causing mitral valve abnormalities with left heart hypoplasia'),

            ('septal_defects', 'left_sided_lesions',
             'Septal & Left Heart Defects',
             'Genes affecting both septal formation and left ventricular development'),

            ('septal_defects', 'aortic_defects',
             'Septal & Aortic Defects',
             'Genes affecting ventricular septum and aortic valve/arch development'),

            ('outflow_tract', 'vascular_defects',
             'Outflow Tract & Vascular Defects',
             'Genes affecting conotruncal and major vessel development'),

            ('outflow_tract', 'left_sided_lesions',
             'Outflow Tract & Left Heart Defects',
             'Genes affecting both outflow tract and left ventricular development'),
        ]

        logger.info(f"Creating intersection groups from {len(intersection_pairs)} biological combinations...")

        # Track which genes have been assigned to intersections
        genes_used_in_intersections = defaultdict(set)  # category -> genes used in intersections
        intersection_count = 0

        # Create intersection groups
        for cat1, cat2, description, rationale in intersection_pairs:
            if cat1 not in categories or cat2 not in categories:
                logger.debug(f"  Skipping {cat1} × {cat2}: category not found")
                continue

            genes1 = categories[cat1]['genes']
            genes2 = categories[cat2]['genes']

            # Find intersection
            intersection_genes = genes1 & genes2

            if len(intersection_genes) > 0:
                sorted_genes = sorted(intersection_genes)
                group_name = f"{cat1}__x__{cat2}"

                final_groups[group_name] = {
                    'description': description,
                    'biological_rationale': rationale,
                    'genes': sorted_genes
                }

                # Track which genes from each category were used
                genes_used_in_intersections[cat1].update(intersection_genes)
                genes_used_in_intersections[cat2].update(intersection_genes)

                intersection_count += 1
                logger.info(f"✓ Intersection {cat1} × {cat2}: {len(intersection_genes)} genes")
            else:
                logger.debug(f"  Skipping {cat1} × {cat2}: no shared genes")

        logger.info(f"Created {intersection_count} intersection groups")

        # Add base categories for remaining genes (not used in any intersection)
        logger.info(f"Creating base groups for genes not in intersections...")

        base_groups_added = 0
        skipped_empty = 0
        skipped_all_used = 0

        for cat_name, cat_info in categories.items():
            genes = cat_info['genes']

            # Skip empty categories
            if len(genes) == 0:
                skipped_empty += 1
                continue

            # Get genes NOT used in any intersection for this category
            remaining_genes = genes - genes_used_in_intersections.get(cat_name, set())

            if len(remaining_genes) == 0:
                logger.debug(f"  Skipping {cat_name}: all genes used in intersections")
                skipped_all_used += 1
                continue

            # Add base group with remaining genes (NO SPLITTING)
            final_groups[cat_name] = {
                'description': cat_info['description'],
                'biological_rationale': cat_info.get('biological_rationale', ''),
                'genes': sorted(remaining_genes)
            }

            if len(remaining_genes) > max_genes_per_group:
                logger.warning(f"⚠ Base group {cat_name}: {len(remaining_genes)} genes (exceeds target of {max_genes_per_group})")
            else:
                logger.info(f"✓ Base group {cat_name}: {len(remaining_genes)} genes")

            base_groups_added += 1

        logger.info(f"Added {base_groups_added} base groups, skipped {skipped_empty + skipped_all_used} categories")

        # Remove groups exceeding max_genes_per_group
        logger.info(f"Removing groups exceeding {max_genes_per_group} genes...")
        groups_before_filter = len(final_groups)

        groups_to_remove = []
        for name, info in final_groups.items():
            if len(info['genes']) > max_genes_per_group:
                groups_to_remove.append((name, len(info['genes'])))

        if groups_to_remove:
            logger.info(f"Removing {len(groups_to_remove)} groups that exceed {max_genes_per_group} genes:")
            for name, size in sorted(groups_to_remove, key=lambda x: x[1], reverse=True):
                logger.info(f"  ✗ Removing {name}: {size} genes")
                del final_groups[name]
        else:
            logger.info(f"✓ All groups within {max_genes_per_group} gene limit")

        # Summary statistics
        total_genes_in_final_groups = sum(len(g['genes']) for g in final_groups.values())
        logger.info(f"Final summary:")
        logger.info(f"  Groups before filtering: {groups_before_filter}")
        logger.info(f"  Groups removed: {len(groups_to_remove)}")
        logger.info(f"  Final groups: {len(final_groups)}")
        logger.info(f"  Total gene assignments: {total_genes_in_final_groups}")

        return final_groups

    def write_gene_groups(
        self,
        groups: Dict[str, Dict],
        output_dir: str,
        general_chd: Set[str],
        max_genes_per_group: int
    ) -> None:
        """
        Write gene groups to CSV files with standardized format.

        Output Format:
        --------------
        Each CSV file contains:
        - Single column named 'x' (for compatibility with downstream tools)
        - One gene symbol per row
        - No index column
        - Genes sorted alphabetically

        File Naming:
        ------------
        Pattern: chd_genes_{category_name}.csv
        Examples:
            - chd_genes_septal_defects.csv
            - chd_genes_outflow_tract_1.csv (for split groups)
            - chd_genes_general.csv

        Parameters:
        -----------
        groups : Dict[str, Dict]
            Final groups with gene lists
        output_dir : str
            Output directory path
        general_chd : Set[str]
            Genes with general CHD annotation
        max_genes_per_group : int
            Maximum genes per file (for splitting general_chd if needed)
        """
        logger.info(f"Writing gene groups to {output_dir}...")

        output_path = Path(output_dir)
        if not output_path.exists():
            logger.info(f"Creating output directory: {output_path}")
            output_path.mkdir(parents=True, exist_ok=True)
        else:
            logger.debug(f"Output directory already exists: {output_path}")

        files_created = 0
        total_genes_written = 0

        # Write specific category groups
        logger.info(f"Writing {len(groups)} category group files...")
        for group_name, group_info in groups.items():
            filename = f"chd_genes_{group_name}.csv"
            filepath = output_path / filename

            # Create single-column dataframe
            df_out = pd.DataFrame({'x': group_info['genes']})
            df_out.to_csv(filepath, index=False)

            files_created += 1
            total_genes_written += len(group_info['genes'])
            logger.info(f"[{files_created:2d}] Created {filename}: "
                       f"{len(group_info['genes'])} genes - {group_info['description']}")
            logger.debug(f"     First gene: {group_info['genes'][0]}, "
                        f"Last gene: {group_info['genes'][-1]}")

        # Write general CHD group (NO SPLITTING, FILTER IF TOO LARGE)
        if general_chd:
            general_genes = sorted(general_chd)
            logger.info(f"Processing general CHD group ({len(general_genes)} genes)...")

            if len(general_genes) > max_genes_per_group:
                logger.warning(f"✗ Skipping general CHD group: {len(general_genes)} genes "
                             f"(exceeds limit of {max_genes_per_group})")
            else:
                filename = "chd_genes_general.csv"
                filepath = output_path / filename
                df_out = pd.DataFrame({'x': general_genes})
                df_out.to_csv(filepath, index=False)

                files_created += 1
                total_genes_written += len(general_genes)
                logger.info(f"[{files_created:2d}] Created {filename}: "
                           f"{len(general_genes)} genes - General CHD")
                logger.debug(f"     First gene: {general_genes[0]}, "
                            f"Last gene: {general_genes[-1]}")

        logger.info(f"✓ Successfully created {files_created} gene group files")
        logger.info(f"  Total genes written: {total_genes_written}")


def main():
    """
    Main execution function for CHD gene grouping pipeline.

    Pipeline Steps:
    ---------------
    1. Parse command-line arguments
    2. Load and filter input data (Human species only)
    3. Categorize genes into biological groups
    4. Split large groups to meet size constraints
    5. Write output CSV files
    6. Generate summary statistics

    Exit Codes:
    -----------
    0: Success
    1: Error (file not found, invalid data, etc.)
    """
    # Configure argument parser with detailed help
    parser = argparse.ArgumentParser(
        description='Create biologically meaningful CHD gene subgroups from scRNA-seq data',
        epilog="""
        Example usage:
          %(prog)s --in data/scRNA_CHD/subtype_spec_genes.csv --max-genes 50
          %(prog)s --in input.csv --max-groups 15 --output-dir results/

        Output:
          Creates multiple CSV files in output directory, each containing
          a group of genes with shared CHD associations. Files are named
          'chd_genes_{category}.csv' with single column 'x' containing
          gene symbols.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--in',
        dest='input_file',
        required=True,
        help='Input CSV file with columns: Gene, Species, Associated CHD'
    )
    parser.add_argument(
        '--max-groups',
        type=int,
        default=20,
        help='Maximum number of groups to create (default: 20, informational only)'
    )
    parser.add_argument(
        '--max-genes',
        type=int,
        default=50,
        help='Maximum genes per group; larger groups will be split (default: 50)'
    )
    parser.add_argument(
        '--output-dir',
        default='data/scRNA_CHD',
        help='Output directory for gene group CSV files (default: data/scRNA_CHD)'
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

    # Log pipeline configuration
    logger.info("=" * 70)
    logger.info("CHD Gene Grouping Pipeline - Starting")
    logger.info("=" * 70)
    logger.info(f"Configuration:")
    logger.info(f"  Input file:       {args.input_file}")
    logger.info(f"  Output directory: {args.output_dir}")
    logger.info(f"  Max groups:       {args.max_groups}")
    logger.info(f"  Max genes/group:  {args.max_genes}")
    logger.info(f"  Verbose mode:     {args.verbose}")
    logger.info("=" * 70)

    # Validate inputs
    logger.info("Step 1: Validating inputs...")
    input_path = Path(args.input_file)
    if not input_path.exists():
        logger.error(f"✗ Input file not found: {args.input_file}")
        sys.exit(1)
    logger.info(f"✓ Input file exists: {args.input_file}")

    if args.max_genes < 1:
        logger.error(f"✗ Invalid --max-genes value: {args.max_genes} (must be ≥ 1)")
        sys.exit(1)
    logger.info(f"✓ Parameters validated")

    try:
        # Initialize grouper
        logger.info("\nStep 2: Initializing CHD gene grouper...")
        grouper = CHDGeneGrouper()
        logger.info("✓ Grouper initialized")

        # Load input data
        logger.info("\nStep 3: Loading input data...")
        logger.info(f"Reading CSV file: {args.input_file}")
        df = pd.read_csv(args.input_file)
        logger.info(f"✓ Loaded {len(df)} total entries")
        logger.debug(f"  Columns: {', '.join(df.columns.tolist())}")

        # Validate required columns
        logger.info("Validating required columns...")
        required_cols = ['Gene', 'Species', 'Associated CHD']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"✗ Missing required columns: {missing_cols}")
            sys.exit(1)
        logger.info(f"✓ All required columns present")

        # Filter for Human species only
        logger.info("\nStep 4: Filtering for Human species...")
        species_counts = df['Species'].value_counts()
        logger.debug(f"Species distribution: {dict(species_counts)}")

        df_human = df[df['Species'] == 'Human'].copy()
        n_genes = df_human['Gene'].nunique()
        logger.info(f"✓ Filtered to {len(df_human)} Human entries with {n_genes} unique genes")
        logger.debug(f"  Excluded {len(df) - len(df_human)} non-human entries")

        if len(df_human) == 0:
            logger.error("✗ No Human entries found in input data")
            sys.exit(1)

        # Categorize genes into biological groups
        logger.info("\nStep 5: Categorizing genes into biological groups...")
        categories, general_chd, gene_assignments = grouper.categorize_genes(
            df_human,
            max_groups=args.max_groups,
            max_genes_per_group=args.max_genes
        )
        logger.info("✓ Gene categorization complete")

        # Create intersection-based groups
        logger.info("\nStep 6: Creating intersection-based groups...")
        final_groups = grouper.create_intersection_groups(
            categories,
            args.max_genes,
            gene_assignments
        )
        logger.info("✓ Intersection-based grouping complete")

        # Check if we exceeded max_groups
        logger.info("\nStep 7: Validating group count...")
        total_groups = len(final_groups) + (1 if general_chd else 0)
        if total_groups > args.max_groups:
            logger.warning(f"⚠ Created {total_groups} groups, exceeding max of {args.max_groups}")
            logger.warning(f"  Consider increasing --max-genes parameter to reduce group count")
        else:
            logger.info(f"✓ Created {total_groups} groups (within max of {args.max_groups})")

        # Write output files
        logger.info("\nStep 8: Writing output files...")
        grouper.write_gene_groups(
            final_groups,
            args.output_dir,
            general_chd,
            args.max_genes
        )
        logger.info("✓ All output files written")

        # Calculate and display summary statistics
        logger.info("\nStep 9: Generating summary statistics...")
        total_genes_assigned = sum(len(g['genes']) for g in final_groups.values()) + len(general_chd)

        # Calculate multi-assignment statistics
        multi_assigned_genes = [gene for gene, cats in gene_assignments.items() if len(cats) > 1]
        logger.debug(f"Genes assigned to multiple categories: {len(multi_assigned_genes)}")
        if multi_assigned_genes and args.verbose:
            logger.debug(f"Examples: {', '.join(multi_assigned_genes[:5])}")

        print(f"\n{'='*78}")
        print(f"{'CHD Gene Grouping Summary':^78}")
        print(f"{'='*78}")
        print(f"  Input file:              {args.input_file}")
        print(f"  Total input entries:     {len(df_human)}")
        print(f"  Unique genes (Human):    {n_genes}")
        print(f"  Groups created:          {total_groups}")
        print(f"  Max genes per group:     {args.max_genes}")
        print(f"  Total gene assignments:  {total_genes_assigned}")
        print(f"  Multi-assigned genes:    {len(multi_assigned_genes)}")
        print(f"  Output directory:        {args.output_dir}")
        print(f"{'='*78}\n")

        logger.info("=" * 70)
        logger.info("✓ Pipeline completed successfully")
        logger.info("=" * 70)

    except Exception as e:
        logger.error(f"Pipeline failed with error: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()


"""
USAGE EXAMPLES
==============

Basic usage:
-----------
python make_chd_gene_groups.py --in data/scRNA_CHD/subtype_spec_genes.csv

With custom parameters:
----------------------
python make_chd_gene_groups.py \\
    --in data/scRNA_CHD/subtype_spec_genes.csv \\
    --max-genes 100 \\
    --output-dir results/chd_groups/

With verbose logging:
--------------------
python make_chd_gene_groups.py \\
    --in data/scRNA_CHD/subtype_spec_genes.csv \\
    --max-genes 50 \\
    --verbose

Output files will be created in the format:
    data/scRNA_CHD/chd_genes_septal_defects_1.csv
    data/scRNA_CHD/chd_genes_outflow_tract_1.csv
    data/scRNA_CHD/chd_genes_pulmonary_defects_1.csv
    ...

Each file contains a single column 'x' with gene symbols:
    x
    GENE1
    GENE2
    GENE3
    ...

BIOLOGICAL CATEGORIES
====================

The script organizes genes into 12 biologically meaningful categories:

1. septal_defects          - Ventricular and atrial septal defects
2. outflow_tract           - Conotruncal and outflow tract defects
3. pulmonary_defects       - Pulmonary valve and artery defects
4. aortic_defects          - Aortic valve and arch defects
5. left_sided_lesions      - Left heart obstructive lesions
6. single_ventricle        - Single ventricle and complex defects
7. vascular_defects        - Patent ductus and vascular anomalies
8. valve_abnormalities     - Atrioventricular valve malformations
9. pulmonary_venous        - Pulmonary venous connection defects
10. laterality_defects     - Heterotaxy and laterality defects
11. right_sided_lesions    - Right ventricular and tricuspid defects
12. coronary_anomalies     - Coronary artery anomalies
13. general                - Non-specific CHD associations

LOGGING LEVELS
==============

Standard logging (default):
- Pipeline step progress (9 steps total)
- Category sizes and assignments
- File creation with gene counts
- Summary statistics

Verbose logging (--verbose):
- Individual gene-CHD type assignments
- Progress percentage during processing
- Species distribution statistics
- Gene ranges for each output file
- Multi-assignment details

EXIT CODES
==========
0 - Success
1 - Error (file not found, missing columns, invalid parameters, etc.)
"""
