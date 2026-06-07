from Bio.SeqUtils.ProtParam import ProteinAnalysis
def charge_isp(seq):
    pa = ProteinAnalysis(seq)
    intestinal_pH = 6.5
    diffeernce = abs(pa.isoelectric_point() - intestinal_pH)
    return [abs(pa.charge_at_pH(intestinal_pH)), diffeernce ]