import pinetree.pinetree as pt

sim = pt.Simulation(120, 1, cell_volume=1.1e-15)

phage = pt.Genome("phage", 200)

phage.add_promoter("p1", 1, 35, {"rnapol": 1e7, "ecolipol": 1e7})
phage.add_terminator("t1", 199, 200, {"rnapol": 1.0, "ecolipol": 1.0})

phage.add_gene("my_gene", 100, 150, 70, 100, 1e7)

phage.add_mask(50, ["rnapol", "ecolipol"])

sim.register_genome(phage)

sim.add_polymerase("rnapol", 35, 50, 10)
sim.add_polymerase("ribosome", 30, 30, 1000)

pt.seed(34)

sim.run("test2")
