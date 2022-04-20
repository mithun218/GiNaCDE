---
title: 'GiNaCDE: the high-performance F-expansion and First Integral Methods with C++ library for solving Nonlinear Differential Equations'
tags:
  - C++
  - symbolic computations
  - nonlinear partial differential equations
  - F-expansion method  
  - First integral method 
authors:
  - name: Mithun Bairagi
    orcid: 0000-0002-9678-4625
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Physics, The University of Burdwan, Golapbag 713104, West Bengal, India
   index: 1
date: 5 January 2021
bibliography: paper.bib

---

# Summary

GiNaCDE is a free and open-source C++ library that solves NLPDEs (Nonlinear Partial Differential Equations) and NLODEs (Nonlinear Ordinary Differential Equations) based on the high-performance F-expansion [@fexpn024], modified F-expansion (mF-expansion) [@modfexpn] and First Integral Methods (FIM) [@fim0]. It unifies and gathers various research on the F-expansion method and FIM [@fexpn024; @fexpn024_1; @fexpn123; @0246; @234; @modfexpn; @fim0; @mirza; @complexTwt1]. It implements two versions of the F-expansion method: F-expansion and modified F-expansion methods. It can be used to get exact traveling-wave solutions of a wide variety of NLPDEs arising in different scientific fields by applying any one of the three available methods (F-expansion, mF-expansion, and first integral method). GiNaCDE is built on a pure C++ symbolic library GiNaC [@ginac], and its algebraic manipulations are performed by the classes provided by the GiNaC library. It also has a rich Graphical User Interface (GUI) which makes it more user-friendly.

The library has been designed to solve the NLPDEs, which have the following general form
\begin{equation}\label{geneq}
    F\left( {\alpha_i,u,u_t,u_{x_1},u_{x_1} \ldots ,u_{x_m},u_{tt},u_{t{x_1}},u_{t{x_2}}, \ldots ,u_{t{x_m}},u_{{x_1}{x_1}},u_{{x_1}{x_2}} \ldots ,u_{{x_1}{x_m}} \ldots } \right) = 0,
\end{equation}
where $t,x_1,x_2, \ldots x_m$ are independent variables, $u=u(t,x_1,x_2, \ldots x_m)$ is a dependent variable, $\alpha_i(i=1,2,\ldots,n)$ are the parameters. Here $F$ must be a polynomial in $u$ and its derivatives. The primary intention of GiNaCDE is to facilitate the development and validation of new exact analytical solutions of NLPDEs and NLODEs. Using the GiNaCDE library, we have successfully solved various kinds of NLPDEs, including higher
nonlinearity terms, higher-derivative terms, and complex NLPDE. Some of them are: Eckhaus equation, seventh-order Sawada-Kotara equations, fifth-order Generalized Korteweg–De Vries (KdV) equation, perturbed nonlinear Schrödinger (NLS) equation with Kerr Law Nonlinearity, and Kudryashov-Sinelshchikov equation. 

# Statement of need 

The NLODEs and NLPDEs play an important role in the theoretical sciences to explain many nonlinear phenomena in various fields of science, such as biology, chemistry, engineering, solid-state physics, plasma physics, optical fibers. The exact (closed-form) traveling-wave solutions of such NLPDEs give much extra information, which helps us to study the result more deeply. The knowledge of closed-form solutions of NLODEs and NLPDEs helps to test the degree of accuracy of numerical solvers and also facilitates stability analysis. In the past few decades, many powerful methods have been presented to seek exact solutions of NLPDEs, such as F-expansion method and first integral method. 
The F-expansion method was first proposed by Zhou et al. [@fexpn024]. The first integral method was first introduced by Feng [@fim0] in solving the Compound Burgers-KdV Equation, which is based on the ring theory of commutative algebra. Later, these methods have been further improved in a number of research works. Following some research works [@fexpn024; @fexpn024_1; @fexpn123; @0246; @234; @modfexpn; @fim0; @mirza; @complexTwt1], GiNaCDE intends to gather and unify the possible combinations of the different revised versions of these methods for improvements of solution procedures. In this context, on the basis of the aforementioned research works of these methods, we have presented the high-performance algorithms of F-expansion, modified F-expansion, and first integral methods in the documentation[^1] file. The GiNaCDE software uses these algorithms to solve the NLPDEs and NLODEs.

In order to solve the NLPDEs, many computer packages are available.
In 1996, Parkes and Duffy [@atfm] had implemented tanh-expansion in their Mathematica package ATFM. Later complete implementation of tanh-expansion has been done by Li and Liu (2002) [@rath] designing the Maple package RATH. Baldwin et al. [@pdespclpkg1] have developed the Mathematica package *PDESpecialSolutions.m* which admits polynomial solutions in tanh, sech, combinations thereof, JacobiSN, JacobiCN. RAEEM [@raeem] is one of the most popular packages written in the Maple programming language, which is a comprehensive and complete implementation of some powerful methods such as the tanh-method, the extended tanh-method, the Jacobi elliptic function method, and the elliptic equation method. One can note that most packages have been developed within commercially available software frameworks, such as Maple and Mathematica. Besides this, all these computer packages have implemented the function-expansion methods. One serious drawback of the function-expansion method is that the solutions which contain functions other than some specific type of functions, such as tanh, sech, JacobiSN, JacobiCN, are not obtained. Additionally, the non-polynomial forms of these particular functions are not obtained also. On the other hand, we have observed that F-expansion, mF-expansion, and FIM are different kinds of methods, which can overcome the limitations of the function-expansion method. To the best of our knowledge, the computer packages implementing F-expansion and first integral methods are not available so far.
Keeping in mind all the above points of view, we have been motivated to develop a free and open-source computer package or C++ library called GiNaCDE, which implements the F-expansion and first integral methods. 

The symbolic manipulations of GiNaCDE depend only on GiNaC [@ginac]. There are several advantages to use GiNaC over other CAS. GiNaC is a free and open-source pure C++ library. It can accept C++ programming language, a general-purpose object-oriented programming (OOP) language, and it is fast like commercially available computer algebra systems.
Besides the library version of GiNaCDE, we have also developed a GUI version of GiNaCDE called GiNaCDE GUI, which facilitates users to solve NLPDEs automatically without writing programming and compilation each time. This GUI version guides us in each step to obtain the output results.

[^1]: https://github.com/mithun218/GiNaCDE/blob/master/doc/documentation.pdf

# References
