import hts

var v:VCF
assert open(v, "../../data/raw/kaviar.vcf")

var afs = new_seq[float32](5) # size doesn't matter. this will be re-sized as needed
var ds = new_string_of_cap(200) # size doesn't matter. this will be re-sized as needed

for rec in v:
  echo rec, " qual:", rec.QUAL, " filter:", rec.FILTER
  var info = rec.info
  # accessing stuff from the INFO field is meant to be as fast as possible, allowing
  # the user to re-use memory as needed.
  #info.strings("DS", ds)
  #info.floats("AF", afs)
  #echo acs 