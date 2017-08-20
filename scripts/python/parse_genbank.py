#! usr/bin/env python3

'''
Parse the T7 genbank record into a pinetree parameter file. The range between
fast and slow codons is configurable via the command line. Whether or not gene 
10 is deoptimized is also configurable via command line arguments.

Many parameter values are hardcoded into this script. I haven't come up with a 
reasonable alternative that doesn't include some sort of meta-parameter file. 
Introducing a meta-parameter file will make managing runs much more confusing, 
so I'm avoiding it. 
'''

from Bio import Entrez, SeqIO
from yaml import dump, safe_load
import math

# CONFIGURATION OPTIONS
RUNTIME = 1800
TIME_STEP = 5
SEED = 34
PHI10_BIND = 1.82e7  # Binding constant for phi10
# Genes to ignore from genbank
IGNORE_GENES = ["gene 10B",
                "possible gene 5.5-5.7",
                "gene 4.1",
                "gene 4B",
                "gene 0.6A",
                "gene 0.6B",
                "possible gene 0.6B",
                "gene 0.5",
                "gene 0.4"]
# Regulatory elements to ignore from genbank
IGNORE_REGULATORY = ["E. coli promoter E[6]",
                     "T7 promoter phiOR",
                     "T7 promoter phiOL",
                     "E. coli promoter A0 (leftward)"]

RELABEL_GENES = {"gene 2": "gp-2",
                 "gene 1": "rnapol-1",
                 "gene 3.5": "lysozyme-3.5",
                 "gene 0.7": "protein_kinase-0.7"}

# Optimal E. coli codons
OPT_CODONS_E_COLI = {'A':['GCT'],
                     'R':['CGT', 'CGC'],
                     'N':['AAC'],
                     'D':['GAC'],
                     'C':['TGC'],
                     'Q':['CAG'],
                     'E':['GAA'],
                     'G':['GGT', 'GGC'],
                     'H':['CAC'],
                     'I':['ATC'],
                     'L':['CTG'],
                     'F':['TTC'],
                     'P':['CCG'],
                     'S':['TCT', 'TCC'],
                     'T':['ACT', 'ACC'],
                     'Y':['TAC'],
                     'V':['GTT', 'GTA']}

def set_up_simulation():
    ''' Set up simulation parameters. '''
    sim = {'seed': SEED,
           'runtime': RUNTIME,
           'time_step': 5,
           'cell_volume': 1.1e-15,
           'debug': False}
    return sim

def set_up_genome():
    ''' Set up genome parameters. '''
    genome = {'name': 'T7',
              'copy_number': 1,
              'entered': 500,
              'mask_interactions':
                  ['rnapol-1',
                   'rnapol-3.5',
                   'ecolipol',
                   'ecolipol-p',
                   'ecolipol-2',
                   'ecolipol-2-p']}

def set_up_polymerases():
    ''' Set up polymerases. '''
    polymerases = [{'copy_number': 0,
                    'footprint': 35,
                    'name': 'rnapol-1',
                    'speed': 230},
                   {'copy_number': 0,
                    'footprint': 35,
                    'name': 'rnapol-3.5',
                    'speed': 230},
                   {'copy_number': 0,
                    'footprint': 35,
                    'name': 'ecolipol',
                    'speed': 45},
                   {'copy_number': 0,
                    'footprint': 35,
                    'name': 'ecolipol-p',
                    'speed': 45},
                   {'copy_number': 0,
                    'footprint': 35,
                    'name': 'ecolipol-2',
                    'speed': 45},
                   {'copy_number': 0,
                    'footprint': 35,
                    'name': 'ecolipol-2-p',
                    'speed': 45}]
    return polymerases

def set_up_ribosomes():
    ''' Set up ribosomes. '''
    ribosomes = [{'binding_constant': '1e7',
                  'copy_number': 0,
                  'footprint': 30,
                  'name': 'ribosome',
                  'speed': 30}]
    return ribosomes

def set_up_species():
    ''' Set up species. '''
    # Assume 10000 total ribosomes and 1800 ecoli polymerases
    species = [{'copy_number': 10000, 'name': 'bound_ribosome'},
               {'copy_number': 1800, 'name': 'bound_ecolipol'},
               {'copy_number': 0, 'name': 'bound_ecolipol_p'},
               {'copy_number': 0, 'name': 'ecoli_genome'},
               {'copy_number': 0, 'name': 'ecoli_transcript'}]
    return species

def set_up_reactions():
    reactions = safe_load('''
    - name: ecoli_transcripts # Assumes that each ribosome binds to one RBS
      propensity: 1e6
      reactants:
          - ecoli_transcript
          - ribosome
      products:
          - bound_ribosome
    - name: ecoli_transcripts_reverse # Assumes that ~1000nt transcript and 30bp/s
      propensity: 0.04
      reactants:
          - bound_ribosome
      products:
          - ribosome
          - ecoli_transcript
    - name: ecoli_transcript_degradation
      propensity: 0.001925  # 1st-order decay based on 6min half-life
      reactants:
          - ecoli_transcript
      products:
          - degraded_transcript
    - name: ecoli_promoters
      propensity: 1e7 # Corresponds to weak E. coli promoter
      reactants:
          - ecolipol
          - ecoli_genome
      products:
          - bound_ecolipol
    - name: ecoli_promoters_p
      propensity: 0.3e7 # Corresponds to weak E. coli promoter
      reactants:
          - ecolipol-p
          - ecoli_genome
      products:
          - bound_ecolipol_p
    - name: ecoli_promoters_reverse # Assumes 1 new RBS per 1000nt gene
      propensity: 0.04
      reactants:
          - bound_ecolipol
      products:
          - ecolipol
          - ecoli_genome
          - ecoli_transcript
    - name: ecoli_promoters_p_reverse # Assumes 1 new RBS per 1000nt gene
      propensity: 0.04
      reactants:
          - bound_ecolipol_p
      products:
          - ecolipol-p
          - ecoli_genome
          - ecoli_transcript
    - name: reaction1
      propensity: 3.8e7
      reactants:
          - protein_kinase-0.7
          - ecolipol
      products:
          - ecolipol-p
          - protein_kinase-0.7
    - name: reaction2
      propensity: 3.8e7
      reactants:
          - protein_kinase-0.7
          - ecolipol-2
      products:
          - ecolipol-2-p
          - protein_kinase-0.7
    - name: reaction3
      propensity: 3.8e7
      reactants:
          - gp-2
          - ecolipol
      products:
          - ecolipol-2
    - name: reaction4
      propensity: 3.8e7
      reactants:
          - gp-2
          - ecolipol-p
      products:
          - ecolipol-2-p
    - name: reaction5
      propensity: 1.1
      reactants:
          - ecolipol-2
      products:
          - gp-2
          - ecolipol
    - name: reaction6
      propensity: 1.1
      reactants:
          - ecolipol-2-p
      products:
          - gp-2
          - ecolipol-p
    - name: reaction7
      propensity: 3.8e7
      reactants:
          - lysozyme-3.5
          - rnapol-1
      products:
          - rnapol-3.5
    - name: reaction8
      propensity: 3.5
      reactants:
          - rnapol-3.5
      products:
          - rnapol-1
          - lysozyme-3.5
    ''')
    return reactions

def set_up():
    data = """
    simulation:
        seed: 34
        runtime: 1800
        time_step: 5
        cell_volume: 1.1e-15
        debug: False
    # Genome parameters
    genome:
        name: T7
        copy_number: 1
        entered: 500
        mask_interactions:
            - rnapol-1
            - rnapol-3.5
            - ecolipol
            - ecolipol-p
            - ecolipol-2
            - ecolipol-2-p

    polymerases:
    - name: rnapol-1
      copy_number: 0
      speed: 230
      footprint: 35
    - name: rnapol-3.5
      copy_number: 0
      speed: 230
      footprint: 35
    - name: ecolipol
      copy_number: 0
      speed: 45
      footprint: 35
    - name: ecolipol-p
      copy_number: 0
      speed: 45
      footprint: 35
    - name: ecolipol-2
      copy_number: 0
      speed: 45
      footprint: 35
    - name: ecolipol-2-p
      copy_number: 0
      speed: 45
      footprint: 35

    ribosomes:
    - name: ribosome
      copy_number: 0
      speed: 30
      footprint: 30
      binding_constant: 1e7

    species:
    - name: bound_ribosome
      copy_number: 10000 # Assume 10000 total ribosomes bound
    - name: bound_ecolipol
      copy_number: 1800
    - name: bound_ecolipol_p
      copy_number: 0
    - name: ecoli_genome
      copy_number: 0
    - name: ecoli_transcript
      copy_number: 0

    reactions:
    - name: ecoli_transcripts # Assumes that each ribosome binds to one RBS
      propensity: 1e6
      reactants:
          - ecoli_transcript
          - ribosome
      products:
          - bound_ribosome
    - name: ecoli_transcripts_reverse # Assumes that ~1000nt transcript and 30bp/s
      propensity: 0.04
      reactants:
          - bound_ribosome
      products:
          - ribosome
          - ecoli_transcript
    - name: ecoli_transcript_degradation
      propensity: 0.001925  # 1st-order decay based on 6min half-life
      reactants:
          - ecoli_transcript
      products:
          - degraded_transcript
    - name: ecoli_promoters
      propensity: 1e7 # Corresponds to weak E. coli promoter
      reactants:
          - ecolipol
          - ecoli_genome
      products:
          - bound_ecolipol
    - name: ecoli_promoters_p
      propensity: 0.3e7 # Corresponds to weak E. coli promoter
      reactants:
          - ecolipol-p
          - ecoli_genome
      products:
          - bound_ecolipol_p
    - name: ecoli_promoters_reverse # Assumes 1 new RBS per 1000nt gene
      propensity: 0.04
      reactants:
          - bound_ecolipol
      products:
          - ecolipol
          - ecoli_genome
          - ecoli_transcript
    - name: ecoli_promoters_p_reverse # Assumes 1 new RBS per 1000nt gene
      propensity: 0.04
      reactants:
          - bound_ecolipol_p
      products:
          - ecolipol-p
          - ecoli_genome
          - ecoli_transcript
    - name: reaction1
      propensity: 3.8e7
      reactants:
          - protein_kinase-0.7
          - ecolipol
      products:
          - ecolipol-p
          - protein_kinase-0.7
    - name: reaction2
      propensity: 3.8e7
      reactants:
          - protein_kinase-0.7
          - ecolipol-2
      products:
          - ecolipol-2-p
          - protein_kinase-0.7
    - name: reaction3
      propensity: 3.8e7
      reactants:
          - gp-2
          - ecolipol
      products:
          - ecolipol-2
    - name: reaction4
      propensity: 3.8e7
      reactants:
          - gp-2
          - ecolipol-p
      products:
          - ecolipol-2-p
    - name: reaction5
      propensity: 1.1
      reactants:
          - ecolipol-2
      products:
          - gp-2
          - ecolipol
    - name: reaction6
      propensity: 1.1
      reactants:
          - ecolipol-2-p
      products:
          - gp-2
          - ecolipol-p
    - name: reaction7
      propensity: 3.8e7
      reactants:
          - lysozyme-3.5
          - rnapol-1
      products:
          - rnapol-3.5
    - name: reaction8
      propensity: 3.5
      reactants:
          - rnapol-3.5
      products:
          - rnapol-1
          - lysozyme-3.5
    """

    return safe_load(data)

def get_promoter_interactions(name):
    '''
    Calculate promoter binding strengths. The relative strengths defined here
    come from 2012 Covert, et al paper.
    '''
    ecoli_strong = ["E. coli promoter A1",
                    "E. coli promoter A2",
                    "E. coli promoter A3"]
    ecoli_weak = ["E. coli B promoter",
                  "E. coli C promoter"]
    phi1_3 = ["T7 promoter phi1.1A",
              "T7 promoter phi1.1B",
              "T7 promoter phi1.3",
              "T7 promoter phi1.5",
              "T7 promoter phi1.6"]
    phi3_8 = ["T7 promoter phi2.5",
              "T7 promoter phi3.8",
              "T7 promoter phi4c",
              "T7 promoter phi4.3",
              "T7 promoter phi4.7"]
    phi6_5 = ["T7 promoter phi6.5"]
    phi9 = ["T7 promoter phi9"]
    phi10 = ["T7 promoter phi10"]
    phi13 = ["T7 promoter phi13",
             "T7 promoter phi17"]

    if name in ecoli_strong:
        return {'ecolipol': {'binding_constant': 10e4},
                'ecolipol-p': {'binding_constant': 3e4}}
    elif name in ecoli_weak:
        return {'ecolipol': {'binding_constant': 1e4},
                'ecolipol-p': {'binding_constant': 0.3e4}}
    elif name in phi1_3:
        return {'rnapol-1': {'binding_constant': PHI10_BIND*0.01},
                'rnapol-3.5': {'binding_constant': PHI10_BIND*0.01*0.5}}
    elif name in phi3_8:
        return {'rnapol-1': {'binding_constant': PHI10_BIND*0.01},
                'rnapol-3.5': {'binding_constant': PHI10_BIND*0.01*0.5}}
    elif name in phi6_5:
        return {'rnapol-1': {'binding_constant': PHI10_BIND*0.05},
                'rnapol-3.5': {'binding_constant': PHI10_BIND*0.05*0.5}}
    elif name in phi9:
        return {'rnapol-1': {'binding_constant': PHI10_BIND*0.2},
                'rnapol-3.5': {'binding_constant': PHI10_BIND*0.2*0.5}}
    elif name in phi10:
        return {'rnapol-1': {'binding_constant': PHI10_BIND},
                'rnapol-3.5': {'binding_constant': PHI10_BIND*0.5}}
    elif name in phi13:
        return {'rnapol-1': {'binding_constant': PHI10_BIND*0.1},
                'rnapol-3.5': {'binding_constant': PHI10_BIND*0.1*0.5}}
    else:
        raise ValueError("Promoter strength for {0} not assigned.".format(name))

def get_terminator_interactions(name):
    '''
    Get terminator efficiencies.
    '''
    if name == "E. coli transcription terminator TE":
        return {'ecolipol': {'efficiency': 1.0},
                'ecolipol-p': {'efficiency': 1.0},
                'rnapol-1': {'efficiency': 0.0},
                'rnapol-3.5': {'efficiency': 0.0}}
    elif name == "T7 transcription terminator Tphi":
        return {'rnapol-1': {'efficiency': 0.85},
                'rnapol-3.5': {'efficiency': 0.85}}
    else:
        return {'name': {'efficiency': 0.0}}

def construct_promoter(feature):
    # Convert to inclusive genomic coordinates
    start = feature.location.start.position + 1
    stop = feature.location.end.position
    name = feature.qualifiers["note"][0]
    if stop - start < 35:
        promoter_start = start - 35
    interactions = get_promoter_interactions(name)
    output = {
        "type": "promoter",
        "name": name,
        "start": promoter_start,
        "stop": stop,
        "interactions": interactions
    }
    return output


def construct_terminator(feature):
    # Convert to inclusive genomic coordinates
    start = feature.location.start.position + 1
    stop = feature.location.end.position
    name = feature.qualifiers["note"][0]
    interactions = get_terminator_interactions(name)
    output = {"type": "terminator",
              "name": name,
              "start": start,
              "stop": stop,
              "interactions": interactions}
    return output

def construct_transcript(feature):
    # Conver to genomic coordinates
    start = feature.location.start.position + 1
    stop = feature.location.end.position
    name = feature.qualifiers["note"][0]
    if name in RELABEL_GENES:
        name = RELABEL_GENES[name]
    # Construct CDS parameters for this gene
    transcript = {"type": "transcript",
                  "name": name,
                  "start": start,
                  "stop": stop,
                  "rbs": -30}
    return transcript

def normalize_weights(weights):
    # Average over all CDSs, which will have non-zero weights
    non_zero = sum(1 if i != 0 else 0 for i in weights)
    mean_weight = sum(weights)/non_zero
    norm_weights = [i/mean_weight for i in weights]
    # Replace non-CDS weights with 1
    norm_weights = [1 if i == 0 else i for i in norm_weights]
    return norm_weights

def compute_cds_weights(record, feature, factor, weights):
    # Grab the gene name
    nuc_seq = feature.location.extract(record).seq
    aa_seq = feature.qualifiers["translation"][0]
    for index, nuc in enumerate(nuc_seq):
        aa_index = int(index / 3)
        codon_start = aa_index * 3
        codon = nuc_seq[codon_start:codon_start + 3]
        genome_index = feature.location.start + index
        if aa_index < len(aa_seq):
            if aa_seq[aa_index] in OPT_CODONS_E_COLI:
                if codon in OPT_CODONS_E_COLI[aa_seq[aa_index]]:
                    weights[genome_index] = factor
                else:
                    weights[genome_index] = 1
        # if feature.qualifiers["protein_id"][0] == 'NP_041998.1':
        #     translation_scale_factors[genome_index] = 0.15
    return weights

def recode_gene(protein_id, weights, record):
    min_weight = min(weights)
    for feature in record.features:
        if feature.type == "CDS":
            if feature.qualifiers["protein_id"][0] == protein_id:
                start = feature.location.start
                stop = feature.location.end
                i = start
                while i < stop:
                    weights[i] = min_weight
                    i += 1
    return weights


def main():
    # NOTE: Read genbank record notes carefully to interpret genomic coordinates
    #
    # "Promoters" are really transcription start sites!!
    # Likewise, "terminators" are really the final position in a transcript
    #
    # All start and stop positions coming form biopython follow PYTHON 
    # conventions. To convert them back to genbank positions, add 1 to the start
    # position.
    # 

    # Download T7 wild-type genbank records
    Entrez.email = "benjamin.r.jack@gmail.com"
    handle = Entrez.efetch(db="nuccore",
                        id=["NC_001604"],
                        rettype="gb",
                        retmode="text")

    records = SeqIO.parse(handle, "genbank")

    output = set_up()

    output["elements"] =[]

    for record in records:

        weights = [0.0]*len(record.seq)

        for feature in record.features:
            # Convert to inclusive genomic coordinates
            start = feature.location.start.position + 1
            stop = feature.location.end.position
            name = ''
            if "note" in feature.qualifiers:
                name = feature.qualifiers["note"][0]
            # Grab promoters and terminators
            if feature.type == "regulatory":
                if name in IGNORE_REGULATORY:
                    continue
                # Construct promoter
                if "promoter" in feature.qualifiers["regulatory_class"]:
                    output["elements"].append(construct_promoter(feature))
                # Construct terminator params
                if "terminator" in feature.qualifiers["regulatory_class"]:
                    output["elements"].append(construct_terminator(feature))
            # Grab genes/CDSes
            if feature.type == "gene":
                if name in IGNORE_GENES:
                    continue
                output['elements'].append(construct_transcript(feature))
            if feature.type == "CDS":
                # Grab the gene name
                weights = compute_cds_weights(record, feature, 10, weights)

    norm_weights = normalize_weights(weights)
    rec_norm_weights = recode_gene("NP_041998.1", norm_weights.copy(), record)
    output["genome"]["translation_weights"] = rec_norm_weights
    output["genome"]["length"] = len(record.seq)

if __name__ == "__main__":
    main()
