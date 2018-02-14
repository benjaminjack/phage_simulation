import pinetree.pinetree as pt

sim = pt.Simulation(run_time=1000, time_step=5, cell_volume=1e-15)

phage = pt.Genome(name, length=39000)

phage.add_promoter(name, start, stop, interactions)

phage.add_terminator(name, start, stop, efficiency)

phage.add_gene(name, start, stop, rbs_start, rbs_stop, rbs_strength)

phage.set_translation_weights(weights=[])

phage.add_mask(entered=500, interactions)

sim.register_genome(phage)

sim.add_ribosome(footprint, mean_speed, copy_number)

sim.add_polymerase(name, footprint, mean_speed, copy_number)

sim.add_species()

sim.add_reaction()

sim.run()
