# ISMRM2022 VSASL Bloch Simulations

This is the code used for the Velocity Selective ASL Bloch simulation tutorial at ISMRM 2022

`run_VSASL_Simulations.m` is the top level file for generating the VSASL modules and running the simulations. 

Bloch simulations are performed using an adapted version of Brian Hargreaves' MEX Bloch simulation code. The original version may be found at http://www-mrsrl.stanford.edu/~brian/blochsim/.

This tutorial was designed and written by [Joseph G. Woods](https://www.ndcn.ox.ac.uk/team/joseph-woods) (University of Oxford), [Dapeng Liu](https://www.researchgate.net/profile/Dapeng-Liu-12) (Johns Hopkins University), [Qin Qin](https://www.hopkinsmedicine.org/profiles/details/qin-qin) (Johns Hopkins University), and  [Yuriko Suzuki](https://www.win.ox.ac.uk/people/yuriko-suzuki) (University of Oxford).

If you use this tutorial for your own teaching, or use some part of the code for any reason, the authors would appreciate that appropriate credit is given.

## Code Contributors

[Joseph G. Woods](https://www.ndcn.ox.ac.uk/team/joseph-woods) (University of Oxford)    
[Dapeng Liu](https://www.researchgate.net/profile/Dapeng-Liu-12) (Johns Hopkins University)    
[Eric C. Wong](https://profiles.ucsd.edu/eric.wong) (University of California San Diego)    
[Jia Guo](https://profiles.ucr.edu/app/home/profile/jiag) (University of California Riverside)    

## References

1. Woods JG, Wong EC, Boyd EC, Bolar DS. VESPA ASL: VElocity and SPAtially Selective Arterial Spin Labeling. Magnetic Resonance in Med. 2022;87(6):2667-2684. doi:https://doi.org/10.1002/mrm.29159
2. Liu D, Li W, Xu F, Zhu D, Shin T, Qin Q. Ensuring both velocity and spatial responses robust to field inhomogeneities for velocity-selective arterial spin labeling through dynamic phase-cycling. Magnetic Resonance in Medicine. 2021;85(5):2723-2734. doi:https://doi.org/10.1002/mrm.28622
3. Holmes JH, Jen ML, Eisenmenger LB, Schubert T, Turski PA, Johnson KM. Spatial dependency and the role of local susceptibility for velocity selective arterial spin labeling (VS-ASL) relative tagging efficiency using accelerated 3D radial sampling with a BIR-8 preparation. Magnetic Resonance in Medicine. 2021;n/a(n/a). doi:https://doi.org/10.1002/mrm.28726
4. Guo J, Das S, Hernandez‐Garcia L. Comparison of velocity‐selective arterial spin labeling schemes. Magnetic Resonance in Medicine. 2021;85(4):2027-2039. doi:https://doi.org/10.1002/mrm.28572
5. Liu D, Xu F, Li W, Zijl PC van, Lin DD, Qin Q. Improved velocity-selective-inversion arterial spin labeling for cerebral blood flow mapping with 3D acquisition. Magnetic Resonance in Medicine. 2020;84(5):2512-2522. doi:https://doi.org/10.1002/mrm.28310
6. Shin T, Qin Q, Park JY, Crawford RS, Rajagopalan S. Identification and reduction of image artifacts in non-contrast-enhanced velocity-selective peripheral angiography at 3T. Magnetic Resonance in Medicine. 2016;76(2):466-477. doi:https://doi.org/10.1002/mrm.25870
7. Qin Q, van Zijl PCM. Velocity-selective-inversion prepared arterial spin labeling. Magnetic Resonance in Medicine. 2016;76(4):1136-1148. doi:https://doi.org/10.1002/mrm.26010
8. Qin Q, Shin T, Schär M, Guo H, Chen H, Qiao Y. Velocity-selective magnetization-prepared non-contrast-enhanced cerebral MR angiography at 3 Tesla: Improved immunity to B0/B1 inhomogeneity. Magnetic Resonance in Medicine. 2016;75(3):1232-1241. doi:https://doi.org/10.1002/mrm.25764
9. Guo J, Meakin J a., Jezzard P, Wong EC. An optimized design to reduce eddy current sensitivity in velocity-selective arterial spin labeling using symmetric BIR-8 pulses. Magnetic Resonance in Medicine. 2015;73(3):1085-1094. doi:https://doi.org/10.1002/mrm.25227
10. Shin T, Worters PW, Hu BS, Nishimura DG. Non-contrast-enhanced renal and abdominal MR angiography using velocity-selective inversion preparation. Magnetic Resonance in Medicine. 2013;69(5):1268-1275. doi:https://doi.org/10.1002/mrm.24356
11. Shin T, Hu BS, Nishimura DG. Off-resonance-robust velocity-selective magnetization preparation for non-contrast-enhanced peripheral MR angiography. Magnetic Resonance in Medicine. 2013;70(5):1229-1240. doi:https://doi.org/10.1002/mrm.24561
12. Meakin JA, Jezzard P. An optimized velocity selective arterial spin labeling module with reduced eddy current sensitivity for improved perfusion quantification. Magnetic Resonance in Medicine. 2013;69(3):832-838. doi:https://doi.org/10.1002/mrm.24302
13. Bolar DS, Rosen BR, Sorensen AG, Adalsteinsson E. QUantitative Imaging of eXtraction of oxygen and TIssue consumption (QUIXOTIC) using venular-targeted velocity-selective spin labeling. Magnetic Resonance in Medicine. 2011;66(6):1550-1562. doi:https://doi.org/10.1002/mrm.22946
14. Wong EC, Guo J. BIR-4 based B1 and B0 insensitive velocity selective pulse trains. In: Proceedings of the 18th Annual Meeting of the ISMRM, Stockholm, Sweden, 2010. 2853. https://cds.ismrm.org/protected/10MProceedings/PDFfiles/2853_7262.pdf.
15. Wu WC, Wong EC. Intravascular effect in velocity-selective arterial spin labeling: The choice of inflow time and cutoff velocity. NeuroImage. 2006;32(1):122-128. doi:https://doi.org/10.1016/j.neuroimage.2006.03.001
16. Wong EC, Cronin M, Wu WC, Inglis B, Frank LR, Liu TT. Velocity-selective arterial spin labeling. Magnetic Resonance in Medicine. 2006;55(6):1334-1341. doi:https://doi.org/10.1002/mrm.20906
17. Norris DG, Schwarzbauer C. Velocity Selective Radiofrequency Pulse Trains. Journal of Magnetic Resonance. 1999;137(1):231-236. doi:https://doi.org/10.1006/jmre.1998.1690



