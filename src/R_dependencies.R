
args = commandArgs(trailingOnly=TRUE)
package = args[1]
print(package)
install.packages(package, repo='http://cran.stat.auckland.ac.nz')
