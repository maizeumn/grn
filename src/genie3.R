#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                     help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                     dest="verbose", help="Print little output")
parser$add_argument("fi", nargs=1, help="Input (expression matrix) file")
parser$add_argument("fo", nargs=1, help="Output file")
parser$add_argument("-p", "--thread", type="integer", default=1,
                     help="Num. threads to use [default %(default)s]")
args <- parser$parse_args()

if( args$verbose ) { 
    write("writing some verbose output to standard error...\n", stderr()) 
}

if( file.access(fi) == -1 ) {
    stop(sprintf("Input file ( %s ) cannot be accessed", fi))
}
