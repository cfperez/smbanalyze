from __future__ import with_statement

def main():
    """
    Reads the .energy and .ct file output of UNAFold and turns it into a more
    legible format for other Python analysis scripts to read.
    """
    import os
    from folding import NAStruct, NAList
    from optparse import OptionParser

    usage = "usage: %prog [options] energy-file [ct-file]"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", action="store_true", dest="verbose")

    (options, args) = parser.parse_args()

    if not args:
	raise ValueError, "Must include energy filename as first positional argument"

    ct_file = ''

    try:
	ct_file = args[1]
    except:
	ct_file = args[0].rpartition('.')[0] + '.ct'
    finally:
	if not os.path.isfile(ct_file):
	    raise IOError, "Couldn't locate ct file: " + ct_file

    with open(ct_file) as ct:
	NA_list = NAList(ct)

    with open(args[0]) as input:
	NA_list.energize(input)

    print NA_list.view()

if __name__ == "__main__":
    main()
