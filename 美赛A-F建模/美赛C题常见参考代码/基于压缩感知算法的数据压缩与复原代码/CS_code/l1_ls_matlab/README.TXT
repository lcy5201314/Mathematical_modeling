l1_ls_matlab.zip 中的程序解决一类最优化问题，使用说明和应用举例都在《l1_ls_useruide》中详细描述。

l1_ls: A Matlab solver for l1-regularized least squares problems.

BETA version, May 10 2008

COPYRIGHT (c) 2007  Kwangmoo Koh, Seung-Jean Kim and Stephen Boyd.

Permission to use, copy, modify, and distribute this software for
any purpose without fee is hereby granted, provided that this entire
notice is included in all copies of any software which is or includes
a copy or modification of this software and in all copies of the
supporting documentation for such software.
This software is being provided "as is", without any express or
implied warranty.  In particular, the authors do not make any
representation or warranty of any kind concerning the merchantability
of this software or its fitness for any particular purpose.

--------------------------------------------------------------------------------

Please see the comments in the source code, for more information.

--------------------------------------------------------------------------------
Files of interest:


SOLVERS
    （主要的函数都在这里）
l1_ls.m : main solver
    
l1_ls_nonneg.m : main solver (with non-negativity constraints)



EXAMPLES
  (自带例子，可以学会这么使用)  
simple_example.m : example for the basic usage
    
operator_example.m : example that illustrates the object programming feature of Matlab



UTILITY
    
find_lambdamax_l1_ls : utility file for finding maximum value of lambda
    find_lambdamax_l1_ls_nonneg : finding maximum value of lambda of
 l1-regularized LSP with nonnegativity constraints
