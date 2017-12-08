from Bio import Entrez, SeqIO
import pinetree.pinetree as pt

CELL_VOLUME = 1.1e-15
PHI10_BIND = 1.82e7  # Binding constant for phi10
MASK_START = 500
IGNORE_REGULATORY = ["E. coli promoter E[6]",
                     "T7 promoter phiOR",
                     "T7 promoter phiOL",
                     "E. coli promoter A0 (leftward)"]

IGNORE_GENES = ["gene 10B",
                "possible gene 5.5-5.7",
                "gene 4.1",
                "gene 4B",
                "gene 0.6A",
                "gene 0.6B",
                "possible gene 0.6B",
                "gene 0.5",
                "gene 0.4"]

RELABEL_GENES = {"gene 2": "gp-2",
                 "gene 1": "rnapol-1",
                 "gene 3.5": "lysozyme-3.5",
                 "gene 0.7": "protein_kinase-0.7"}

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
        return {'ecolipol': 10e4,
                'ecolipol-p': 3e4}
    elif name in ecoli_weak:
        return {'ecolipol': 1e4,
                'ecolipol-p': 0.3e4}
    elif name in phi1_3:
        return {'rnapol-1': PHI10_BIND * 0.01,
                'rnapol-3.5': PHI10_BIND * 0.01 * 0.5}
    elif name in phi3_8:
        return {'rnapol-1': PHI10_BIND * 0.01,
                'rnapol-3.5': PHI10_BIND * 0.01 * 0.5}
    elif name in phi6_5:
        return {'rnapol-1': PHI10_BIND * 0.05,
                'rnapol-3.5': PHI10_BIND * 0.05}
    elif name in phi9:
        return {'rnapol-1': PHI10_BIND * 0.2,
                'rnapol-3.5': PHI10_BIND * 0.2}
    elif name in phi10:
        return {'rnapol-1': PHI10_BIND,
                'rnapol-3.5': PHI10_BIND}
    elif name in phi13:
        return {'rnapol-1': PHI10_BIND * 0.1,
                'rnapol-3.5': PHI10_BIND * 0.1}
    else:
        raise ValueError(
            "Promoter strength for {0} not assigned.".format(name))


def construct_promoter(feature):
    # Convert to inclusive genomic coordinates
    start = feature.location.start.position + 1
    stop = feature.location.end.position
    name = feature.qualifiers["note"][0]
    length = stop - start
    if length < 35:
        start = start - 35
    tracker = pt.SpeciesTracker.get_instance()
    if stop < MASK_START:
        tracker.increment(name, 1)
    interactions = get_promoter_interactions(name)
    prom = pt.Promoter(name, start, stop, list(interactions.keys()))
    bind_reactions = list()
    for partner, constant in interactions.items():
        if partner == 'rnapol-1' or partner == 'rnapol-3.5':
            pol = pt.Polymerase(partner, 35, 230)
        else:
            pol = pt.Polymerase(partner, 35, 45)
        bind = pt.Bind(constant, CELL_VOLUME, name, pol)
        bind_reactions.append(bind)
    return prom, bind_reactions


def get_terminator_interactions(name):
    '''
    Get terminator efficiencies.
    '''
    if name == "E. coli transcription terminator TE":
        return {'ecolipol': 1.0,
                'ecolipol-p': 1.0,
                'rnapol-1': 0.0,
                'rnapol-3.5': 0.0}
    elif name == "T7 transcription terminator Tphi":
        return {'rnapol-1': 0.85,
                'rnapol-3.5': 0.85}
    else:
        return {'name': 0.0}

def construct_terminator(feature):
    # Convert to inclusive genomic coordinates
    start = feature.location.start.position + 1
    stop = feature.location.end.position
    name = feature.qualifiers["note"][0]
    interactions = get_terminator_interactions(name)
    term = pt.Terminator(name, start, stop, list(interactions.keys()), interactions)
    return term

def construct_transcript(feature):
    # Conver to genomic coordinates
    start = feature.location.start.position + 1
    stop = feature.location.end.position
    name = feature.qualifiers["note"][0]
    if name in RELABEL_GENES:
        name = RELABEL_GENES[name]
    # Construct CDS parameters for this gene
    rbs = pt.Promoter("rbs", start - 30, start, ["ribosome"])
    rbs.gene = name
    stop_site = pt.Terminator("tstop", stop - 1, stop, 
                           ["ribosome"], {"ribosome": 1.0})
    stop_site.reading_frame = start % 3
    stop_site.gene = name
    return rbs, stop_site

def main():
    sim = pt.Simulation()

    sim.stop_time = 1500  # ~22 minutes
    sim.time_step = 5

    # Download T7 wild-type genbank records
    Entrez.email = "benjamin.r.jack@gmail.com"
    handle = Entrez.efetch(db="nuccore",
                           id=["NC_001604"],
                           rettype="gb",
                           retmode="text")

    records = SeqIO.parse(handle, "genbank")

    dna_elements = list()
    transcript_template = list()

    for record in records:
        genome_length = len(record.seq)
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
                    prom, bind_reactions = construct_promoter(feature)
                    for bind in bind_reactions:
                        sim.register_reaction(bind)
                    dna_elements.append(prom)
                # Construct terminator params
                if "terminator" in feature.qualifiers["regulatory_class"]:
                    dna_elements.append(construct_terminator(feature))
            # Grab genes/CDSes
            if feature.type == "gene":
                if name in IGNORE_GENES:
                    continue
                rbs, stop_site = construct_transcript(feature)
                transcript_template.append(rbs)
                transcript_template.append(stop_site)
    
    mask = pt.Mask("mask", MASK_START, genome_length, ["rnapol-1",
                                                "rnapol-3.5",
                                                "ecolipol",
                                                "ecolipol-p",
                                                "ecolipol-2",
                                                "ecolipol-2-p"])
    
    genome = pt.Genome("phage", genome_length, dna_elements, 
                       transcript_template, mask)
    
    sim.register_genome(genome)

    tracker = pt.SpeciesTracker.get_instance()
    tracker.increment("rnapol-1", 0)
    tracker.increment("rnapol-3.5", 0)
    tracker.increment("ecolipol", 0)
    tracker.increment("ecolipol-p", 0)
    tracker.increment("ecolipol-2", 0)
    tracker.increment("ecolipol-2-p", 0)

    ribosome = pt.Polymerase("ribosome", 30, 30)
    ribo_bind = pt.Bind(1e7, CELL_VOLUME, "rbs", ribosome)
    sim.register_reaction(ribo_bind)
    tracker.increment("ribosome", 0)
    tracker.increment("bound_ribosome", 10000)

    tracker.increment("bound_ecolipol", 1800)
    tracker.increment("bound_ecolipol_p", 0)
    tracker.increment("ecoli_genome", 0)
    tracker.increment("ecoli_transcript", 0)

    sim.register_reaction(pt.SpeciesReaction(1e6,
                                            CELL_VOLUME,
                                            ["ecoli_transcript", "ribosome"],
                                            ["bound_ribosome"]))
    
    sim.register_reaction(pt.SpeciesReaction(0.04,
                                            CELL_VOLUME,
                                            ["bound_ribosome"],
                                            ["ribosome", "ecoli_transcript"]))

    sim.register_reaction(pt.SpeciesReaction(0.001925,
                                            CELL_VOLUME,
                                            ["ecoli_transcript"],
                                            ["degraded_transcript"]))
    
    sim.register_reaction(pt.SpeciesReaction(1e7,
                                            CELL_VOLUME,
                                            ["ecolipol", "ecoli_genome"],
                                            ["bound_ecolipol"]))

    sim.register_reaction(pt.SpeciesReaction(0.3e7,
                                            CELL_VOLUME,
                                            ["ecolipol-p", "ecoli_genome"],
                                            ["bound_ecolipol_p"]))
    
    sim.register_reaction(pt.SpeciesReaction(0.04,
                                            CELL_VOLUME,
                                            ["bound_ecolipol"],
                                            ["ecolipol", "ecoli_genome", "ecoli_transcript"]))
    
    sim.register_reaction(pt.SpeciesReaction(0.04,
                                            CELL_VOLUME,
                                            ["bound_ecolipol_p"],
                                            ["ecolipol-p", "ecoli_genome", "ecoli_transcript"]))
    
    sim.register_reaction(pt.SpeciesReaction(3.8e7,
                                            CELL_VOLUME,
                                            ["protein_kinase-0.7", "ecolipol"],
                                            ["ecolipol-p", "protein_kinase-0.7"]))
    
    sim.register_reaction(pt.SpeciesReaction(3.8e7,
                                            CELL_VOLUME,
                                            ["protein_kinase-0.7", "ecolipol-2"],
                                            ["ecolipol-2-p", "protein_kinase-0.7"]))
    
    sim.register_reaction(pt.SpeciesReaction(3.8e7,
                                            CELL_VOLUME,
                                            ["gp-2", "ecolipol"],
                                            ["ecolipol-2"]))

    sim.register_reaction(pt.SpeciesReaction(3.8e7,
                                            CELL_VOLUME,
                                            ["gp-2", "ecolipol-p"],
                                            ["ecolipol-2-p"]))
    
    sim.register_reaction(pt.SpeciesReaction(1.1,
                                            CELL_VOLUME,
                                            ["ecolipol-2-p"],
                                            ["gp-2", "ecolipol-p"]))
    
    sim.register_reaction(pt.SpeciesReaction(1.1,
                                            CELL_VOLUME,
                                            ["ecolipol-2"],
                                            ["gp-2", "ecolipol"]))
    
    sim.register_reaction(pt.SpeciesReaction(3.8e9,
                                            CELL_VOLUME,
                                            ["lysozyme-3.5", "rnapol-1"],
                                            ["rnapol-3.5"]))

    sim.register_reaction(pt.SpeciesReaction(3.5,
                                            CELL_VOLUME,
                                            ["rnapol-3.5"],
                                            ["lysozyme-3.5", "rnapol-1"]))
    
    pt.seed(34)

    sim.run("test")
    
if __name__ == "__main__":
    main()

