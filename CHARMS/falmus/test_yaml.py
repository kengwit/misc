import yaml


filename = 'Mandible.yaml'
fh = open(filename, 'r')
doc = yaml.load(fh) 
fh.close()


print(doc['algorithms']['gsubmeshes'])