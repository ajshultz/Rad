"""
Copyright (c) 2012, 2013 The PyPedia Project, http://www.pypedia.com
<br>All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: 

# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

http://www.opensource.org/licenses/BSD-2-Clause
"""

__pypdoc__ = """
Method: Pairwise_linkage_disequilibrium
Link: http://www.pypedia.com/index.php/Pairwise_linkage_disequilibrium
Retrieve date: Thu, 21 Nov 2013 01:35:50 +0200



Computes the [http://en.wikipedia.org/wiki/Linkage_disequilibrium linkage disequilibrium] between two SNPs. Computes the r-squared, D' and the estimated haplotype frequencies. It has been developed to produce the same values as plink with --ld option. (http://pngu.mgh.harvard.edu/~purcell/plink/ld.shtml)

'''Parameters:'''
* '''SNP_1''' : A tuple with the alleles of the first SNP. For example (("A","A"), ("A","G"), ("G","G"))
* '''SNP_2''' : A tuple with the alleles of the secons SNP. for example (("A","T"), ("A","A"), ("A","A"))

[[Category:Algorithms]]
[[Category:Bioinformatics]]
[[Category:Validated]]

* See also: http://rannala.org/books/CUPChap3.pdf



"""

def Minor_allele_frequency(
        allele1=None,
        allele2=None,
        observations_chr1 = None,
        observations_chr2 = None):

        """
        >>> Minor_allele_frequency("A", "C", ["A", "A", "C", "C", "C"], ["A", "C", "C", "C","C"])
        (0.3, 'A', 10)
        """

        allele1_chr1 = observations_chr1.count(allele1)
        allele2_chr1 = observations_chr1.count(allele2)

        if observations_chr2:
                allele1_chr2 = observations_chr2.count(allele1)
                allele2_chr2 = observations_chr2.count(allele2)
        else:
                allele1_chr2 = 0
                allele2_chr2 = 0

        observations = allele1_chr1 + allele2_chr1 + allele1_chr2 + allele2_chr2
        allele1_maf = float((allele1_chr1 + allele1_chr2)) / observations

        if allele1_maf < 0.5:
                return (allele1_maf, allele1, observations)
        else:
                return (1.0 - allele1_maf, allele2, observations)


def SNP_alleles(
	SNP = None,
	missing = '0',
):
	allele_1 = None

	for gen in SNP:
		if gen[0] != missing:
			if not allele_1:
				allele_1 = gen[0]
			elif gen[0] != allele_1:
				return (allele_1, gen[0])
		if gen[1] != missing:
			if not allele_1:
				allele_1 = gen[1]
			elif gen[1] != allele_1:
				return (allele_1, gen[1])

	return (allele_1, allele_1)

def Pairwise_linkage_disequilibrium(
	SNP_1 = None,
	SNP_2 = None,
	):
	
	"""
	>>> print Pairwise_linkage_disequilibrium([('A','A'), ('A','G'), ('G','G'), ('G','A')], [('A','A'), ('A','G'), ('G','G'), ('A','A')]) 
	{'haplotypes': [('AA', 0.499999999973935, 0.3125), ('AG', 2.606498430642265e-11, 0.1875), ('GA', 0.12500000002606498, 0.3125), ('GG', 0.374999999973935, 0.1875)], 'R_sq': 0.5999999998331841, 'Dprime': 0.9999999998609868}
	"""

	ret = {}

	(SNP_1_allele_1, SNP_1_allele_2) = SNP_alleles(SNP_1)
	(SNP_2_allele_1, SNP_2_allele_2) = SNP_alleles(SNP_2)

	Haplotype_1_1 = [x[0] for x in SNP_1]
	Haplotype_1_2 = [x[1] for x in SNP_1]
	Haplotype_2_1 = [x[0] for x in SNP_2]
	Haplotype_2_2 = [x[1] for x in SNP_2]

	MAF_1 = Minor_allele_frequency(SNP_1_allele_1, SNP_1_allele_2, Haplotype_1_1, Haplotype_1_2)
	MAF_2 = Minor_allele_frequency(SNP_2_allele_1, SNP_2_allele_2, Haplotype_2_1, Haplotype_2_2)

	if MAF_1[1] == SNP_1_allele_1:
		freq_1_1 = MAF_1[0]
		freq_1_2 = 1.0 - freq_1_1
	elif MAF_1[1] == SNP_1_allele_2:
		freq_1_2 = MAF_1[0]
		freq_1_1 = 1.0 - freq_1_2

	if MAF_2[1] == SNP_2_allele_1:
		freq_2_1 = MAF_2[0]
		freq_2_2 = 1.0 - freq_2_1
	elif MAF_2[1] == SNP_2_allele_2:
		freq_2_2 = MAF_2[0]
		freq_2_1 = 1.0 - freq_2_2

	freq_A_B = freq_1_1 * freq_2_1
	freq_A_b = freq_1_1 * freq_2_2
	freq_a_B = freq_1_2 * freq_2_1
	freq_a_b = freq_1_2 * freq_2_2

	Haps = zip(Haplotype_1_1, Haplotype_1_2, Haplotype_2_1, Haplotype_2_2)


	N11 = len([1 for h11, h12, h21, h22 in Haps if   h11 == SNP_1_allele_1 and h12 == SNP_1_allele_1                                                        and   h21 == SNP_2_allele_1 and h22 == SNP_2_allele_1])
	N21 = len([1 for h11, h12, h21, h22 in Haps if ((h11 == SNP_1_allele_1 and h12 == SNP_1_allele_2) or (h11 == SNP_1_allele_2 and h12 == SNP_1_allele_1)) and   h21 == SNP_2_allele_1 and h22 == SNP_2_allele_1])
	N31 = len([1 for h11, h12, h21, h22 in Haps if   h11 == SNP_1_allele_2 and h12 == SNP_1_allele_2                                                        and   h21 == SNP_2_allele_1 and h22 == SNP_2_allele_1])
	N12 = len([1 for h11, h12, h21, h22 in Haps if   h11 == SNP_1_allele_1 and h12 == SNP_1_allele_1                                                        and ((h21 == SNP_2_allele_1 and h22 == SNP_2_allele_2) or (h21 == SNP_2_allele_2 and h22 == SNP_2_allele_1))])
	N22 = len([1 for h11, h12, h21, h22 in Haps if ((h11 == SNP_1_allele_1 and h12 == SNP_1_allele_2) or (h11 == SNP_1_allele_2 and h12 == SNP_1_allele_1)) and ((h21 == SNP_2_allele_1 and h22 == SNP_2_allele_2) or (h21 == SNP_2_allele_2 and h22 == SNP_2_allele_1))])
	N32 = len([1 for h11, h12, h21, h22 in Haps if   h11 == SNP_1_allele_2 and h12 == SNP_1_allele_2                                                        and ((h21 == SNP_2_allele_1 and h22 == SNP_2_allele_2) or (h21 == SNP_2_allele_2 and h22 == SNP_2_allele_1))])
	N13 = len([1 for h11, h12, h21, h22 in Haps if   h11 == SNP_1_allele_1 and h12 == SNP_1_allele_1                                                        and   h21 == SNP_2_allele_2 and h22 == SNP_2_allele_2])
	N23 = len([1 for h11, h12, h21, h22 in Haps if ((h11 == SNP_1_allele_1 and h12 == SNP_1_allele_2) or (h11 == SNP_1_allele_2 and h12 == SNP_1_allele_1)) and   h21 == SNP_2_allele_2 and h22 == SNP_2_allele_2])
	N33 = len([1 for h11, h12, h21, h22 in Haps if   h11 == SNP_1_allele_2 and h12 == SNP_1_allele_2                                                        and   h21 == SNP_2_allele_2 and h22 == SNP_2_allele_2])

	N1_ = N11 + N12 + N13
	N2_ = N21 + N22 + N23
	N3_ = N31 + N32 + N33

	N_1 = N11 + N21 + N31
	N_2 = N12 + N22 + N32
	N_3 = N13 + N23 + N33

	N = N1_ + N2_ + N3_

	pAB = 0.25
	pab = 0.25
	pAb = 0.25
	paB = 0.25

	for nn in xrange(10):
		nAB = (2*N11) + N12 + N21 + (( (pAB*pab) / ((pAb*paB) + (pAB*pab)) ) * N22)
		nab = (2*N33) + N23 + N32 + (( (pAB*pab) / ((pAb*paB) + (pAB*pab)) ) * N22)
		nAb = (2*N13) + N12 + N23 + (( (pAb*paB) / ((pAb*paB) + (pAB*pab)) ) * N22)
		naB = (2*N31) + N21 + N32 + (( (pAb*paB) / ((pAb*paB) + (pAB*pab)) ) * N22)

		pAB = nAB / (2.0 * N)
		pab = nab / (2.0 * N)
		pAb = nAb / (2.0 * N)
		paB = naB / (2.0 * N)

	ret["haplotypes"] = [
		(SNP_1_allele_1 + SNP_2_allele_1, pAB, freq_A_B),
		(SNP_1_allele_1 + SNP_2_allele_2, pAb, freq_A_b),
		(SNP_1_allele_2 + SNP_2_allele_1, paB, freq_a_B),
		(SNP_1_allele_2 + SNP_2_allele_2, pab, freq_a_b),
		]

	Observed = [pAB, pAb, paB, pab]
	Expected = [freq_A_B, freq_A_b, freq_a_B, freq_a_b]

#	R_sq = sum([((o-e)*(o-e))/e for o,e in zip(Observed, Expected)])
#	ret["R_sq"] = R_sq
	

	pA = pAB + pAb
	pB = pAB + paB

	D = pAB - (pA*pB)

#Altered R_sq calculations so that they match formula from John Wakeley's Coalescent Theory book, and the output from PLINK.	
	R_sq = (D**2)/(pA*(1-pA)*pB*(1-pB)) 
	ret["R_sq"] = R_sq


	if D>=0:
		Dmax = min(pA * (1-pB), (1-pA)*pB)
	else:
		Dmax = min(pA * pB, (1-pA)*(1-pB))

	Dprime = D / Dmax
	ret["Dprime"] = Dprime

	return ret


#Method name =Pairwise_linkage_disequilibrium()
if __name__ == '__main__':
    print __pypdoc__
    SNP_1 = [('A','A'), ('A','G'), ('G','G'), ('G','A')]
    SNP_2 = [('A','A'), ('A','G'), ('G','G'), ('A','A')]

    returned = Pairwise_linkage_disequilibrium(SNP_1=SNP_1, SNP_2=SNP_2)
    if returned:
        print 'Method returned:'
        print str(returned)
