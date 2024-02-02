
try:
    with open("test_ecoli_spades.gfa", mode="r") as graph, open("gplas_input/test_links.txt", mode="w") as links, open("gplas_input/test_raw_nodes.fasta", mode="w") as nodes, open("gplas_input/test_contigs.fasta", mode="w") as contigs:
        for line in graph:
            if line[0] == "S":
                entry = line.split("\t")
                number = "".join([entry[0], entry[1]])
                sequence = entry[2]
                if len(entry) > 4: #unicycler
                    information = "_".join(entry[3:])
                else: #spades
                    information = entry[3]
                nodes.write(f">{number}_{information}{sequence}\n") # "information" already ends with a "\n"
                if len(sequence) >= 1000:
                    contigs.write(f">{number}_{information}{sequence}\n")
            elif line[0] == "L":
                links.write(line)
except IOError as e:
    print("Operation failed:", e.strerror)