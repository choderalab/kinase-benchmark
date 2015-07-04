

def main(verbose, pH):
    import os
    import commands
    from schrodinger import structure # Requires Schrodinger Suite

    # Write input for epik.
    if verbose:
        print "Converting input file to Maestro format..."
    reader = structure.StructureReader("epik-input.mol2")
    writer = structure.StructureWriter("epik-input.mae")
    for st in reader:
        writer.append(st)
    reader.close()
    writer.close()

    # Run epik to enumerate protomers/tautomers and get associated state penalties.
    if verbose:
        print "Running Epik..."
    cmd = '%s/epik -imae epik-input.mae -omae epik-output.mae -pht 10.0 -ms 100 -nt -pKa_atom -ph %f -WAIT' % (os.environ['SCHRODINGER'], pH)
    output = commands.getoutput(cmd)
    if verbose:
        print output

    # Convert output from epik from .mae to .sdf.
    if verbose:
        print "Converting output file to SDF..."
    reader = structure.StructureReader("epik-output.mae")
    writer = structure.StructureWriter("epik-output.sdf")
    for st in reader:
        writer.append(st)
    reader.close()
    writer.close()

    # Also convert to .mol2.
    if verbose:
        print "Converting output file to MOL2..."
    reader = structure.StructureReader("epik-output.mae")
    writer = structure.StructureWriter("epik-output.mol2")
    for st in reader:
        writer.append(st)
    reader.close()
    writer.close()


if __name__ == '__main__':
    import sys

    verbose = bool(sys.argv[1])
    pH = float(sys.argv[2])

    main(verbose, pH)
