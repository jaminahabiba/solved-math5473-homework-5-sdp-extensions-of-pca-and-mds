Download Link: https://assignmentchef.com/product/solved-math5473-homework-5-sdp-extensions-of-pca-and-mds
<br>
<ol>

 <li><em>RPCA</em>: Construct a random rank-<em>r </em>matrix: let <em>A </em>∈ R<em><sup>m</sup></em><sup>×<em>n </em></sup>with <em>a<sub>ij </sub></em>∼ N(0<em>,</em>1) whose top-<em>r </em>singular value/vector is <em>λ<sub>i</sub></em>, <em>u<sub>i </sub></em>∈R<em><sup>m </sup></em>and <em>v<sub>i </sub></em>∈R<em><sup>n </sup></em>(<em>i </em>= 1<em>,…,r</em>), define. Construct a sparse matrix <em>E </em>with <em>p </em>percentage (<em>p </em>∈ [0<em>,</em>1]) nonzero entries distributed uniformly. Then define</li>

</ol>

<em>M </em>= <em>L </em>+ <em>E.</em>

<ul>

 <li>Set <em>m </em>= <em>n </em>= 20, <em>r </em>= 1, and <em>p </em>= 0<em>.</em>1, use Matlab toolbox CVX to formulate a semidefinite program for Robust PCA of <em>M</em>:</li>

</ul>

(trace(<em>W</em><sub>1</sub>) + trace(<em>W</em><sub>2</sub>)) + <em>λ</em>k<em>S</em>k<sub>1                                                                  </sub>(1)

where you can use the matlab implementation in lecture notes as a reference;

<ul>

 <li>Choose different parameters <em>p </em>∈ [0<em>,</em>1] to explore the probability of successful recover;</li>

 <li>Increase <em>r </em>to explore the probability of successful recover;</li>

 <li><em><sup>? </sup></em>Increase <em>m </em>and <em>n </em>to values beyond 50 will make CVX difficult to solve. In this case, use the Augmented Lagrange Multiplier method, e.g. in E. J. Candes, X. Li, Y. Ma, and J. Wright (2009) “Robust Principal Component Analysis?”. Journal of ACM, 58(1), 1-37. Make a code yourself (just a few lines of Matlab or Python) and test it for <em>m </em>= <em>n </em>= 1000. A convergence criterion often used can be for example).</li>

</ul>

<ol start="2">

 <li><em>SPCA</em>: Define three hidden factors:</li>

</ol>

<em>,</em>

where <em>V</em><sub>1</sub><em>,V</em><sub>2</sub>, and <em> </em>are independent. Construct 10 observed variables as follows

<em>,</em>

with <em>j </em>= 1 for <em>i </em>= 1<em>,…,</em>4, <em>j </em>= 2 for <em>i </em>= 5<em>,…,</em>8, and <em>j </em>= 3 for <em>i </em>= 9<em>,</em>10 and  independent for <em>j </em>= 1<em>,</em>2<em>,</em>3, <em>i </em>= 1<em>,…,</em>10.

The first two principal components should be concentrated on (<em>X</em><sub>1</sub><em>,X</em><sub>2</sub><em>,X</em><sub>3</sub><em>,X</em><sub>4</sub>) and (<em>X</em><sub>5</sub><em>,X</em><sub>6</sub><em>,X</em><sub>7</sub><em>,X</em><sub>8</sub>), respectively. This is an example given by H. Zou, T. Hastie, and R. Tibshirani, Sparse principal component analysis, J. Comput. Graphical Statist., 15 (2006), pp. 265-286.

1

<em>Homework 5. SDP Extensions of PCA and MDS                                                                                                                   </em>2

<ul>

 <li>Compute the true covariance matrix Σ (and the sample covariance matrix with <em>n </em>examples, say <em>n </em>= 1000);</li>

 <li>Compute the top 4 principal components of Σ using eigenvector decomposition (byMatlab or R);</li>

 <li>Use Matlab CVX toolbox to compute the first <em>sparse </em>principal component by solving the SDP problem</li>

</ul>

max       trace(Σ<em>X</em>) − <em>λ</em>k<em>X</em>k<sub>1</sub>

<em>s.t.        </em>trace(<em>X</em>) = 1

Choose <em>λ </em>= 0 and other positive numbers to compare your results with normal PCA;

<ul>

 <li>Remove the first sparse PCA from Σ and compute the second sparse PCA with the samecode;</li>

 <li>Again compute the 3rd and the 4th sparse PCA of Σ and compare them against thenormal PCAs.</li>

 <li><em><sup>? </sup></em>Construct an example with 200 observed variables which is hard to deal with by CVX. In this case, try the Augmented Lagrange Multiplier method by Allen Yang et al. (UC Berkeley) whose Matlab codes can be found at <a href="https://www.eecs.berkeley.edu/~yang/software/SPCA/SPCA_ALM.zip">http://www.eecs.berkeley.edu/ </a><a href="https://www.eecs.berkeley.edu/~yang/software/SPCA/SPCA_ALM.zip"><sub>~</sub></a><a href="https://www.eecs.berkeley.edu/~yang/software/SPCA/SPCA_ALM.zip">yang/software/SPCA/SPCA_ALM.zip</a><a href="https://www.eecs.berkeley.edu/~yang/software/SPCA/SPCA_ALM.zip">,</a> or Python scikit-learn Sparse PCA package.</li>

</ul>

<ol start="3">

 <li><em><sup>?</sup>Protein Folding: </em>Consider the 3D structure reconstruction based on incomplete MDS with uncertainty. Data file: <a href="https://yao-lab.github.io/data/protein3D.zip">http://yao-lab.github.io/data/protein3D.zip</a></li>

</ol>

Figure 1: 3D graphs of file PF00018 2HDA.pdf (YES HUMAN/97-144, PDB 2HDA)

In the file, you will find 3D coordinates for the following three protein families:

PF00013 (PCBP1 HUMAN/281-343, PDB 1WVN),

<em>Homework 5. SDP Extensions of PCA and MDS                                                                                                                   </em>3

PF00018 (YES HUMAN/97-144, PDB 2HDA), and PF00254 (O45418 CAEEL/24-118, PDB 1R9H).

For example, the file PF00018 2HDA.pdb contains the 3D coordinates of alpha-carbons for a particular amino acid sequence in the family, YES HUMAN/97-144, read as

VALYDYEARTTEDLSFKKGERFQIINNTEGDWWEARSIATGKNGYIPS where the first line in the file is

97 V 0.967 18.470 4.342

Here

<ul>

 <li>‘97’: start position 97 in the sequence</li>

 <li>‘V’: first character in the sequence</li>

 <li>[<em>x,y,z</em>]: 3D coordinates in unit ˚<em>A</em>.</li>

</ul>

Figure 1 gives a 3D representation of its structure.

Given the 3D coordinates of the amino acids in the sequence, one can computer pairwise distance between amino acids, [<em>d<sub>ij</sub></em>]<em><sup>l</sup></em><sup>×<em>l </em></sup>where <em>l </em>is the sequence length. A <em>contact map </em>is defined to be a graph <em>G<sub>θ </sub></em>= (<em>V,E</em>) consisting <em>l </em>vertices for amino acids such that and edge (<em>i,j</em>) ∈ <em>E </em>if <em>d<sub>ij </sub></em>≤ <em>θ</em>, where the threshold is typically <em>θ </em>= 5˚<em>A </em>or 8˚<em>A </em>here.

Can you recover the 3D structure of such proteins, up to an Euclidean transformation (rotation and translation), given noisy pairwise distances restricted on the contact map graph <em>G<sub>θ</sub></em>, i.e. given noisy pairwise distances between vertex pairs whose true distances are no more than <em>θ</em>? Design a noise model (e.g. Gaussian or uniformly bounded) for your experiments.

When <em>θ </em>= ∞ without noise, classical MDS will work; but for a finite <em>θ </em>with noisy measurements, SDP approach can be useful. You may try the matlab package SNLSDP by Kim-Chuan Toh, Pratik Biswas, and Yinyu Ye, or the facial reduction speed up by Nathan Krislock and Henry Wolkowicz.