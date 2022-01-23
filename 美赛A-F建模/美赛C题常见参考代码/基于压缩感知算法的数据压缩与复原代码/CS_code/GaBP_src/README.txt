此程序与文章Polynomial Linear Programming with Gaussian Belief Propagation.pdf相对应;
资料来源：Gaussian Belief Propagation Resources.htm

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gaussian belief propagation Matlab package

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Writen by Danny Bickson, HUJI & IBM Haifa Lab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Acknowledgment: the following people contributed code/their help
for this package:
N. Sommer - TAU, and E. N. Hoch, HUJI: helped in implementing the LDLC decoder code
A. Zymnis - Stanford: Wrote the NUM simulation
NBP/LDLC encoding matrices where kindly provided by Marilynn Green, Nokia Siemens Networks Research Technology Platforms Dallas, TX.
S. Joshi - Stanford, Wrote the original sensor selection simulation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The following files are basic implementations of the
GaBP algorithm.

1) gabp.m - GaBP algorithm for dense matrices, parallel version

algorithm described in:
Linear Detection via Belief Propagation
By Danny Bickson, Danny Dolev, Ori Shental, Paul H. Siegel and Jack K. Wolf.
In the 45th Annual Allerton Conference on Communication, Control and Computing, Allerton House, Illinois, 
Sept. 07'
1a) run_gabp.m - example script for running gabp.m
2) asynch_GBP.m - GaBP algorithm for dense matrices, serial (asynch. version)
3) sparse_gabp.m - GaBP algorithm for sparse matrices, optimized.
This code was tested with instances of size 500,000 x 500,000 with 4% non zeros
3a) run_sparse_gabp.m - script file for running sprase_gabp.m
4) gabpms.m - Min-Sum algorithm
C.C. Moallemi and B. Van Roy. 揅onvergence of the min-sum algorithm for convex optimization.?
4a) run_gabpms.m - script for running the gabpms.m file
5) gabp_inv.m - runs GaBP for inverting a matrix
5a) run_gabp_inv.m - script for running gabp_inv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The following directories include samples of problems solved using GaBP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1) LDLC - extended low density lattice decoder
This codes needs the KDE Matlab package by A. T. Ihler, found in http://ttic.uchicago.edu/~ihler/code/
After unpacking, you should addpath() in Matlab to the root kde directory.

2) GaBP_convergence_fix - fixing the convergence of the GaBP algorithm
"Fixing the convergence of the GaBP algorithm" by J. K. Johnson, D. Bickson and D. Dolev
In IEEE Internatioanal Synposium on Information Theory, Seoul, South Korea, July 2009. 
http://arxiv.org/abs/0901.4192

3) LP - linear programming example
"Polynomial Linear Programming with Gaussian Belief Propagation. By Danny Bickson, Yoav Tock, 
Ori Shental, Paul H. Seigel, Jack K. Wolf and Danny Dolev. In the Forty-Sixth Annual Allerton 
Conference on Communication, Control, and Computing, Sept. 2008, Allerton House, Monticello, Illinois.
http://arxiv.org/abs/0810.1631

4) NUM - network utility maximization
"Distributed Large Scale Network Utility Maximization", by D. Bickson, Y. Tock, A. Zymnis, S. Boyd and D. Dolev.
In IEEE Internatioanal Synposium on Information Theory, Seoul, South Korea, July 2009.  
http://arxiv.org/abs/0901.2684

5) NBP - non-parametric belief propagation implementation (using sampling and fft() ) - more efficient.

6) ILP - non-parametric belief propagation implementation (using KDE matlab tookbox)

7) Sensor selection example.
Distributed sensor selection via Gaussian belief propagation. D. Bickson and D. Dolev.
In the 47h Annual Allerton Conference on Communication, Control and Computing, Allerton House, Illinois,
Sept. 2009, submitted for publication. 
http://arxiv.org/abs/0907.0931

8) NBP decoder - "Low density lattice decoder via non-parametric belief propagation"
By D. Bickson, A. T. Ihler, H. Avissar and D. Dolver, Submited for publication.
http://arxiv.org/abs/0901.3197



