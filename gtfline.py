def attr_interpreter(attr):
    attr_list = {}
    attr_s = attr.replace("  "," ").split(" ")
    l = int(len(attr_s)/2)
    for i in range(l):
        attr_list[attr_s[2*i]] = attr_s[2*i+1].strip(";").strip("\"")
    return attr_list

class gtf_line:
    def __init__(self, g_line):
        line_tem = g_line.strip().split("\t")
        self.chrom = line_tem[0]
        self.source = line_tem[1]
        self.type = line_tem[2]
        self.start = int(line_tem[3])
        self.end = int(line_tem[4])
        self.score = line_tem[5]
        self.strand = line_tem[6]
        self.phase = line_tem[7]
        self.attributes = attr_interpreter(line_tem[8])
