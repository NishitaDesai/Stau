# SOFTSUSY SUGRA calculation
# SOFTSUSY3.3.7 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.3.7       # version number
Block MODSEL  # Select model
     1    1   # sugra
Block SMINPUTS             # Standard Model inputs
     1    1.27916000e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.18400000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.18000000e+00   # mb(mb)
     6    1.73500000e+02   # Mtop(pole)
     7    1.77699000e+00   # Mtau(pole)
Block MINPAR               # SUSY breaking input parameters
     3    1.00000000e+01   # tanb
     4    1.00000000e+00   # sign(mu)
     1    1.42000000e+02   # m0
     2    7.20000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.63062868e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=1 Desired accuracy=1.00000000e-03 Achieved accuracy=2.10032129e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03912477e+01   # MW
        25     1.17255502e+02   # h0
        35     9.99792844e+02   # H0
        36     9.99594995e+02   # A0
        37     1.00302826e+03   # H+
   1000021     1.61117678e+03   # ~g
   1000022     2.700e+02   # ~neutralino(1)
   1000023     5.67589148e+02   # ~neutralino(2)
   1000024     5.67844618e+02   # ~chargino(1)
   1000025    -8.77817527e+02   # ~neutralino(3)
   1000035     8.88992634e+02   # ~neutralino(4)
   1000037     8.88863793e+02   # ~chargino(2)
   1000001     1.47508664e+03   # ~d_L
   1000002     1.47311729e+03   # ~u_L
   1000003     1.47508323e+03   # ~s_L
   1000004     1.47311387e+03   # ~c_L
   1000005     1.35406831e+03   # ~b_1
   1000006     1.14378743e+03   # ~t_1
   1000011     5.03643262e+02   # ~e_L
   1000012     4.97238090e+02   # ~nue_L
   1000013     5.03688327e+02   # ~mu_L
   1000014     4.97233683e+02   # ~numu_L
   1000015     2.710e+02   # ~stau_1
   1000016     4.95676678e+02   # ~nu_tau_L
   2000001     1.41030308e+03   # ~d_R
   2000002     1.41600446e+03   # ~u_R
   2000003     1.41029950e+03   # ~s_R
   2000004     1.41600081e+03   # ~c_R
   2000005     1.40496211e+03   # ~b_2
   2000006     1.38874804e+03   # ~t_2
   2000011     3.08200531e+02   # ~e_R
   2000013     3.08186043e+02   # ~mu_R
   2000015     6.03663248e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05444319e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97805756e-01   # N_{1,1}
  1  2    -9.15325200e-03   # N_{1,2}
  1  3     6.01570860e-02   # N_{1,3}
  1  4    -2.60962878e-02   # N_{1,4}
  2  1     2.14059529e-02   # N_{2,1}
  2  2     9.81024853e-01   # N_{2,2}
  2  3    -1.57822543e-01   # N_{2,3}
  2  4     1.10562508e-01   # N_{2,4}
  3  1    -2.36950728e-02   # N_{3,1}
  3  2     3.42107648e-02   # N_{3,2}
  3  3     7.05246795e-01   # N_{3,3}
  3  4     7.07739447e-01   # N_{3,4}
  4  1    -5.80000169e-02   # N_{4,1}
  4  2     1.90620250e-01   # N_{4,2}
  4  3     6.88549292e-01   # N_{4,3}
  4  4    -6.97280281e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.74893431e-01   # U_{1,1}
  1  2    -2.22671953e-01   # U_{1,2}
  2  1     2.22671953e-01   # U_{2,1}
  2  2     9.74893431e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.87644717e-01   # V_{1,1}
  1  2    -1.56709646e-01   # V_{1,2}
  2  1     1.56709646e-01   # V_{2,1}
  2  2     9.87644717e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.44871894e-01   # F_{11}
  1  2     9.38649763e-01   # F_{12}
  2  1     9.38649763e-01   # F_{21}
  2  2    -3.44871894e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.85717547e-01   # F_{11}
  1  2     1.68407001e-01   # F_{12}
  2  1    -1.68407001e-01   # F_{21}
  2  2     9.85717547e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     9.81562366e-02   # F_{11}
  1  2     9.95171017e-01   # F_{12}
  2  1     9.95171017e-01   # F_{21}
  2  2    -9.81562366e-02   # F_{22}
Block gauge Q= 1.22378388e+03  # SM gauge couplings
     1     3.63183547e-01   # g'(Q)MSSM DRbar
     2     6.41182545e-01   # g(Q)MSSM DRbar
     3     1.04730075e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.22378388e+03  
  3  3     8.49221073e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.22378388e+03  
  3  3     1.29368514e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.22378388e+03  
  3  3     1.00268030e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.22378388e+03 # Higgs mixing parameters
     1     8.72248358e+02    # mu(Q)MSSM DRbar
     2     9.63362083e+00    # tan beta(Q)MSSM DRbar
     3     2.43726196e+02    # higgs vev(Q)MSSM DRbar
     4     1.03512800e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.22378388e+03  # MSSM DRbar SUSY breaking parameters
     1     3.06347807e+02      # M_1(Q)
     2     5.63451507e+02      # M_2(Q)
     3     1.56831325e+03      # M_3(Q)
    21     2.15225620e+05      # mH1^2(Q)
    22    -7.30779486e+05      # mH2^2(Q)
    31     4.94446685e+02      # meL(Q)
    32     4.94442260e+02      # mmuL(Q)
    33     4.93102949e+02      # mtauL(Q)
    34     2.99351966e+02      # meR(Q)
    35     2.99337048e+02      # mmuR(Q)
    36     2.94792455e+02      # mtauR(Q)
    41     1.42517595e+03      # mqL1(Q)
    42     1.42517248e+03      # mqL2(Q)
    43     1.31631953e+03      # mqL3(Q)
    44     1.36951546e+03      # muR(Q)
    45     1.36951176e+03      # mcR(Q)
    46     1.13234300e+03      # mtR(Q)
    47     1.36266098e+03      # mdR(Q)
    48     1.36265734e+03      # msR(Q)
    49     1.35625654e+03      # mbR(Q)
Block au Q= 1.22378388e+03  
  1  1    -1.59115008e+03      # Au(Q)MSSM DRbar
  2  2    -1.59114332e+03      # Ac(Q)MSSM DRbar
  3  3    -1.23696306e+03      # At(Q)MSSM DRbar
Block ad Q= 1.22378388e+03  
  1  1    -1.93693767e+03      # Ad(Q)MSSM DRbar
  2  2    -1.93693139e+03      # As(Q)MSSM DRbar
  3  3    -1.81305024e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.22378388e+03  
  1  1    -4.24193402e+02      # Ae(Q)MSSM DRbar
  2  2    -4.24185959e+02      # Amu(Q)MSSM DRbar
  3  3    -4.21932840e+02      # Atau(Q)MSSM DRbar
#         PDG            Width
DECAY   1000015     1.3150000E-18   # stau_1 decays
#          BR         NDA      ID1         ID2       ID3       ID4
     1.00E+00          3     1000022        16      -213   # BR(stau_1 -> chi_1 + neu_tau + rho(770))

