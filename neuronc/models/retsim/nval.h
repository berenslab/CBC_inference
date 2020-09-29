/* nval.h */

/* To add parameters, edit "maknval.cc", compile and run "maknval > nval.n",
   then copy nval.n to "nval.h". Remove the param defs from the end of nval.n.
   Remove the nval.n table at the beginning of nval.h.
   Copy nval.h to "nval_var.h", "nval_var.cc", and "nval_var_set.cc"
   and edit these files, then remove this content from nval.h.
   Last, "make clean" and "make retsim".
*/

#define XCONE        0	/* cones */
#define XROD         1	/* rods */
#define HBAT         2	/* hbat */
#define HA           3	/* Type A horizontal cells */
#define HB           4	/* Type B horizontal cells */
#define RBP          5	/* Rod bipolar cells */
#define DBP1         6	/* Depolarizing cone bipolar cell, type 1 */
#define DBP2         7	/* Depolarizing cone bipolar cell, type 2 */
#define DBP3         8	/* Depolarizing cone bipolar cell, type 3 */
#define DBP4         9	/* Depolarizing cone bipolar cell, type 4 */
#define HBP1        10	/* Hyperpolarizing bipolar cell, type 1 */
#define HBP2        11	/* Hyperpolarizing bipolar cell, type 2 */
#define A17         12	/* A17 amacrine cells, feedback to RBP */
#define AII         13	/* AII amacrine cells */
#define SBAC        14	/* Starburst amacrine cells */
#define AM          15	/* Amacrine cell type 1 */
#define AM2         16	/* Amacrine cell type 2 */
#define AM3         17	/* Amacrine cell type 3 */
#define AM4         18	/* Amacrine cell type 4 */
#define AMH         19	/* Amacrine cells, hyperpolarizing */
#define AMH2        20	/* Amacrine cells, hyperpolarizing */
#define AMS         21	/* Amacrine cells, small-field */
#define AMHS        22	/* Amacrine cells, small-field hyperpol */
#define GCA         23	/* Ganglion cells, On-type, alpha */
#define GCB         24	/* Ganglion cells, On-type, beta */
#define DSGC        25	/* Direction-selective ganglion cells */
#define GCAOFF      26	/* Ganglion cells, Off-type, alpha */
#define GCBOFF      27	/* Ganglion cells, Off-type, beta */
#define NCELTYPES   28	/* Number of cell types */

#define MAKE         0	 /* whether to make this cell type */
#define MAKE_DEND    1	 /* whether to make dendrites */
#define MAKE_AXON    2	 /* whether to make axon */
#define MAKE_DIST    3	 /* whether to make axon distal */
#define NMADE        4	 /* number of cells made */
#define MAXNUM       5	 /* maximum number of cells of this type */
#define NCOLOR       6	 /* color of this cell type for display */
#define MAXCOV       7	 /* max coverage factor (for arrays) */
#define MAXSYNI      8	 /* max number of syn input cells */
#define MAXSYNO      9	 /* max number of syn output cells */
#define DENS        10	 /* density of this type (per mm2) */
#define REGU        11	 /* regularity (mean/stdev) of spacing */
#define MORPH       12	 /* morphology (=0 -> file, or artificial) */
#define COMPLAM     13	 /* compartment size (default=complam) */
#define BIOPHYS     14	 /* add biophys properties (chan dens file) */
#define CHNOISE     15	 /* add membrane channel noise properties   */
#define RATIOK      16	 /* set K density values as ratio from Na */
#define VSTART      17	 /* initial resting potential */
#define VREV        18	 /* membrane potential for Rm (VCl) */
#define NRM         19	 /* the cell's Rm, 0 => use default (drm) */
#define NRI         20	 /* the cell's Ri, 0 => use default (dri) */
#define SOMADIA     21	 /* Soma diameter */
#define SOMAZ       22	 /* Z location (x,y loc determ. by array) */
#define DENDARB     23	 /* type of dendritic tree */
#define DENDARBZ    24	 /* dendritic arborization level */
#define DENZDIST    25	 /* dendritic arborization z tolerance */
#define STRATDIA    26	 /* stratif. annulus dia (fract of treedia) */
#define DTIPDIA     27	 /* diameter of dendritic tips */
#define DTREEDIA    28	 /* diameter of dendritic tree */
#define ARBSCALE    29	 /* scale for dia of real morph dend tree */
#define DENDDIA     30	 /* dend dia scale for real morph */
#define AXARBT      31	 /* type of axonal tree */
#define AXARBZ      32	 /* axonal arborization level */
#define AXTIPDIA    33	 /* diameter of axonal tips */
#define AXARBDIA    34	 /* diameter of axonal arbor */
#define MAXSDIST    35	 /* maximum synaptic distance */
#define TAPERSPC    36	 /* space constant of diameter taper */
#define TAPERABS    37	 /* abs diameter for taper */
#define NDENDR      38	 /* number of first-order dendrites */
#define GROWTHR     39	 /* distance thresh for growth of dendrites */
#define SEGLEN      40	 /* length of dendrite segments */

#define CELPRE1     41	 /* cell type to connect to (neg, no conn) */
#define CONPRE1     42	 /* connection number of presyn cell */
#define CELCONV1    43	 /* number of presyn cells to connect to */
#define GROWPOST1   44	 /* grow when making conn from presyn cell */
#define CELPRE2     45	 /* cell type to connect to (neg, no conn) */
#define CONPRE2     46	 /* connection number of presyn cell */
#define CELCONV2    47	 /* number of presyn cells to connect to */
#define GROWPOST2   48	 /* grow when making conn from presyn cell */
#define CELPRE3     49	 /* cell type to connect to (neg, no conn) */
#define CONPRE3     50	 /* connection number of presyn cell */
#define CELCONV3    51	 /* number of presyn cells to connect to */
#define GROWPOST3   52	 /* grow when making conn from presyn cell */
#define CELPRE4     53	 /* cell type to connect to (neg, no conn) */
#define CONPRE4     54	 /* connection number of presyn cell */
#define CELCONV4    55	 /* number of presyn cells to connect to */
#define GROWPOST4   56	 /* grow when making conn from presyn cell */
#define CELPRE5     57	 /* cell type to connect to (neg, no conn) */
#define CONPRE5     58	 /* connection number of presyn cell */
#define CELCONV5    59	 /* number of presyn cells to connect to */
#define GROWPOST5   60	 /* grow when making conn from presyn cell */
#define CELPRE6     61	 /* cell type to connect to (neg, no conn) */
#define CONPRE6     62	 /* connection number of presyn cell */
#define CELCONV6    63	 /* number of presyn cells to connect to */
#define GROWPOST6   64	 /* grow when making conn from presyn cell */
#define CELPRE7     65	 /* cell type to connect to (neg, no conn) */
#define CONPRE7     66	 /* connection number of presyn cell */
#define CELCONV7    67	 /* number of presyn cells to connect to */
#define GROWPOST7   68	 /* grow when making conn from presyn cell */
#define CELPRE8     69	 /* cell type to connect to (neg, no conn) */
#define CONPRE8     70	 /* connection number of presyn cell */
#define CELCONV8    71	 /* number of presyn cells to connect to */
#define GROWPOST8   72	 /* grow when making conn from presyn cell */
#define CELPRE9     73	 /* cell type to connect to (neg, no conn) */
#define CONPRE9     74	 /* connection number of presyn cell */
#define CELCONV9    75	 /* number of presyn cells to connect to */
#define GROWPOST9   76	 /* grow when making conn from presyn cell */
#define CELPRE10    77	 /* cell type to connect to (neg, no conn) */
#define CONPRE10    78	 /* connection number of presyn cell */
#define CELCONV10   79	 /* number of presyn cells to connect to */
#define GROWPOST10  80	 /* grow when making conn from presyn cell */

#define CELPOST1    81	 /* cell type to connect to (neg, no conn) */
#define CONPOST1    82	 /* connection number for postsyn cell */
#define CELDIV1     83	 /* number of postsyn cells to connect to */
#define GROWPRE1    84	 /* grow when making conn to postsyn cell */
#define SYNREG1     85	 /* synaptic region in presyn dendritic tree */
#define SYNREGP1    86	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC1    87	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI1    88	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO1    89	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI1    90	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO1    91	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG1     92	 /* angle for postsynaptic cell */
#define SYNRNG1     93	 /* range of angles for postsynaptic cell */
#define USEDYAD1    94	 /* synapse is dyad using preexisting type */
#define DYADTYP1    95	 /* type of dyad synapse to connect with */
#define AUTAPSE1    96	 /* synapse back to presynaptic node */
#define SYNNUM1     97	 /* number of synapses per connection */
#define SENSCA1     98	 /* synaptic release sensitivity calcium */
#define SRRPOOL1    99	 /* synaptic readily releasable pool */
#define SRRPOOLG1  100	 /* synaptic readily releasable pool gain */
#define SMRRPOOL1  101	 /* synaptic max readily releasable pool */
#define SMAXRATE1  102	 /* maximum sustained synaptic release rate */
#define SGAIN1     103	 /* synaptic gain */
#define SVGAIN1    104	 /* synaptic vgain */
#define SDURH1     105	 /* synaptic high pass time const. */
#define SNFILTH1   106	 /* synaptic high pass nfilt */
#define SHGAIN1    107	 /* synaptic high pass gain */
#define SHOFFS1    108	 /* synaptic high pass offset */
#define SVSIZ1     109	 /* synaptic vesicle size */
#define SCOND1     110	 /* synaptic conductance */
#define SCMUL1     111	 /* synaptic conductance mult for region */
#define SCGRAD1    112	 /* synaptic conductance gradient from soma */
#define SEGRAD1    113	 /* synaptic conductance expon grad fr soma */
#define STHRESH1   114	 /* synaptic threshold */
#define SVNOISE1   115	 /* 1 -> vesicle noise, override, vnoise=0 */
#define SCOV1      116	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR1      117	 /* synaptic event time const. */
#define SFALL1     118	 /* synaptic event fall time const. */
#define SNFILT1    119	 /* synaptic vesicle nfilt */
#define STRCONC1   120	 /* synaptic transmitter concentration. */
#define SRESP1     121	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA1      122	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX1   123	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM1     124	 /* synaptic postsyn Ca pump Km. */
#define SCNFILT1   125	 /* second mesng. nfilt */
#define SCDUR1     126	 /* second mesng. time const. */
#define SCGAIN1    127	 /* synaptic second messenger gain */
#define SCOFF1     128	 /* synaptic second messenger offset */
#define SCNOISE1   129	 /* 1 -> channel noise, override cnoise=0 */
#define SNCHAN1    130	 /* number of channels */
#define SUNIT1     131	 /* synaptic channel unitary conductace */
#define SVREV1     132	 /* synaptic reversal potential */

#define CELPOST2   133	 /* cell type to connect to (neg, no conn) */
#define CONPOST2   134	 /* connection number for postsyn cell */
#define CELDIV2    135	 /* number of postsyn cells to connect to */
#define GROWPRE2   136	 /* grow when making conn to postsyn cell */
#define SYNREG2    137	 /* synaptic region in presyn dendritic tree */
#define SYNREGP2   138	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC2   139	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI2   140	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO2   141	 /* outer rad of annulus in presyn dend tree */
#define SYNANPI2   142	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO2   143	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG2    144	 /* angle for postsynaptic cell */
#define SYNRNG2    145	 /* range of angles for postsynaptic cell */
#define USEDYAD2   146	 /* synapse is dyad using preexisting type */
#define DYADTYP2   147	 /* type of dyad synapse to connect with */
#define AUTAPSE2   148	 /* synapse back to presynaptic node */
#define SYNNUM2    149	 /* number of synapses per connection */
#define SENSCA2    150	 /* synaptic release sensitivity calcium */
#define SRRPOOL2   151	 /* synaptic readily releasable pool */
#define SRRPOOLG2  152	 /* synaptic readily releasable pool gain */
#define SMRRPOOL2  153	 /* synaptic max readily releasable pool */
#define SMAXRATE2  154	 /* maximum sustained synaptic release rate */
#define SGAIN2     155	 /* synaptic gain */
#define SVGAIN2    156	 /* synaptic vgain */
#define SDURH2     157	 /* synaptic high pass time const. */
#define SNFILTH2   158	 /* synaptic high pass nfilt */
#define SHGAIN2    159	 /* synaptic high pass gain */
#define SHOFFS2    160	 /* synaptic high pass offset */
#define SVSIZ2     161	 /* synaptic vesicle size */
#define SCOND2     162	 /* synaptic conductance */
#define SCMUL2     163	 /* synaptic conductance mult for region */
#define SCGRAD2    164	 /* synaptic conductance grad from soma */
#define SEGRAD2    165	 /* synaptic conductance expon grad fr soma */
#define STHRESH2   166	 /* synaptic threshold */
#define SVNOISE2   167	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV2      168	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR2      169	 /* synaptic event time const. */
#define SFALL2     170	 /* synaptic event fall time const. */
#define SNFILT2    171	 /* synaptic vesicle nfilt */
#define STRCONC2   172	 /* synaptic transmitter concentration. */
#define SRESP2     173	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA2      174	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX2   175	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM2     176	 /* synaptic postsyn Ca pump Km. */
#define SCNFILT2   177	 /* second mesng. nfilt */
#define SCDUR2     178	 /* second mesng. time const. */
#define SCGAIN2    179	 /* synaptic second messenger gain */
#define SCOFF2     180	 /* synaptic second messenger offset */
#define SCNOISE2   181	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN2    182	 /* number of channels */
#define SUNIT2     183	 /* synaptic channel unitary conductace */
#define SVREV2     184	 /* synaptic reversal potential */

#define CELPOST3   185	 /* cell type to connect to (neg, no conn) */
#define CONPOST3   186	 /* connection number for postsyn cell */
#define CELDIV3    187	 /* number of postsyn cells to connect to */
#define GROWPRE3   188	 /* grow when making conn to postsyn cell */
#define SYNREG3    189	 /* synaptic region in presyn dendritic tree */
#define SYNREGP3   190	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC3   191	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI3   192	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO3   193	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI3   194	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO3   195	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG3    196	 /* angle for postsynaptic cell */
#define SYNRNG3    197	 /* range of angles for postsynaptic cell */
#define USEDYAD3   198	 /* synapse is dyad using preexisting type */
#define DYADTYP3   199	 /* type of dyad synapse to connect with */
#define AUTAPSE3   200	 /* synapse back to presynaptic node */
#define SYNNUM3    201	 /* number of synapses per connection */
#define SENSCA3    202	 /* synaptic release sensitivity calcium */
#define SRRPOOL3   203	 /* synaptic readily releasable pool */
#define SRRPOOLG3  204	 /* synaptic readily releasable pool gain */
#define SMRRPOOL3  205	 /* synaptic max readily releasable pool */
#define SMAXRATE3  206	 /* maximum sustained synaptic release rate */
#define SGAIN3     207	 /* synaptic gain */
#define SVGAIN3    208	 /* synaptic vgain */
#define SDURH3     209	 /* synaptic high pass time const. */
#define SNFILTH3   210	 /* synaptic high pass nfilt */
#define SHGAIN3    211	 /* synaptic high pass gain */
#define SHOFFS3    212	 /* synaptic high pass offset */
#define SVSIZ3     213	 /* synaptic vesicle size */
#define SCOND3     214	 /* synaptic conductance */
#define SCMUL3     215	 /* synaptic conductance mult for region */
#define SCGRAD3    216	 /* synaptic conductance gradient from soma */
#define SEGRAD3    217	 /* synaptic conductance expon grad fr soma */
#define STHRESH3   218	 /* synaptic threshold */
#define SVNOISE3   219	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV3      220	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR3      221	 /* synaptic event time const. */
#define SFALL3     222	 /* synaptic event fall time const. */
#define SNFILT3    223	 /* synaptic vesicle nfilt */
#define STRCONC3   224	 /* synaptic transmitter concentration. */
#define SRESP3     225	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA3      226	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX3   227	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM3     228	 /* synaptic postsyn Ca pump Km. */
#define SCNFILT3   229	 /* second mesng. nfilt */
#define SCDUR3     230	 /* second mesng. time const. */
#define SCGAIN3    231	 /* synaptic second messenger gain */
#define SCOFF3     232	 /* synaptic second messenger offset */
#define SCNOISE3   233	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN3    234	 /* number of channels */
#define SUNIT3     235	 /* synaptic channel unitary conductace */
#define SVREV3     236	 /* synaptic reversal potential */

#define CELPOST4   237	 /* cell type to connect to (neg, no conn) */
#define CONPOST4   238	 /* connection number for postsyn cell */
#define CELDIV4    239	 /* number of postsyn cells to connect to */
#define GROWPRE4   240	 /* grow when making conn to postsyn cell */
#define SYNREG4    241	 /* synaptic region in presyn dendritic tree */
#define SYNREGP4   242	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC4   243	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI4   244	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO4   245	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI4   246	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO4   247	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG4    248	 /* angle for postsynaptic cell */
#define SYNRNG4    249	 /* range of angles for postsynaptic cell */
#define USEDYAD4   250	 /* synapse is dyad using preexisting type */
#define DYADTYP4   251	 /* type of dyad synapse to connect with */
#define AUTAPSE4   252	 /* synapse back to presynaptic node */
#define SYNNUM4    253	 /* number of synapses per connection */
#define SENSCA4    254	 /* synaptic release sensitivity calcium */
#define SRRPOOL4   255	 /* synaptic readily releasable pool */
#define SRRPOOLG4  256	 /* synaptic readily releasable pool gain */
#define SMRRPOOL4  257	 /* synaptic max readily releasable pool */
#define SMAXRATE4  258	 /* maximum sustained synaptic release rate */
#define SGAIN4     259	 /* synaptic gain */
#define SVGAIN4    260	 /* synaptic vgain */
#define SDURH4     261	 /* synaptic high pass time const. */
#define SNFILTH4   262	 /* synaptic high pass nfilt */
#define SHGAIN4    263	 /* synaptic high pass gain */
#define SHOFFS4    264	 /* synaptic high pass offset */
#define SVSIZ4     265	 /* synaptic vesicle size */
#define SCOND4     266	 /* synaptic conductance */
#define SCMUL4     267	 /* synaptic conductance mult for region */
#define SCGRAD4    268	 /* synaptic conductance gradient from soma */
#define SEGRAD4    269	 /* synaptic conductance expon grad fr soma */
#define STHRESH4   270	 /* synaptic threshold */
#define SVNOISE4   271	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV4      272	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR4      273	 /* synaptic event time const. */
#define SFALL4     274	 /* synaptic event fall time const. */
#define SNFILT4    275	 /* synaptic vesicle nfilt */
#define STRCONC4   276	 /* synaptic transmitter concentration. */
#define SRESP4     277	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA4      278	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX4   279	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM4     280	 /* synaptic postsyn Ca pump Km. */
#define SCNFILT4   281	 /* second mesng. nfilt */
#define SCDUR4     282	 /* second mesng. time const. */
#define SCGAIN4    283	 /* synaptic second messenger gain */
#define SCOFF4     284	 /* synaptic second messenger offset */
#define SCNOISE4   285	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN4    286	 /* number of channels */
#define SUNIT4     287	 /* synaptic channel unitary conductace */
#define SVREV4     288	 /* synaptic reversal potential */

#define CELPOST5   289	 /* cell type to connect to (neg, no conn) */
#define CONPOST5   290	 /* connection number for postsyn cell */
#define CELDIV5    291	 /* number of postsyn cells to connect to */
#define GROWPRE5   292	 /* grow when making conn to postsyn cell */
#define SYNREG5    293	 /* synaptic region in presyn dendritic tree */
#define SYNREGP5   294	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC5   295	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI5   296	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO5   297	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI5   298	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO5   299	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG5    300	 /* angle for postsynaptic cell */
#define SYNRNG5    301	 /* range of angles for postsynaptic cell */
#define USEDYAD5   302	 /* synapse is dyad using preexisting type */
#define DYADTYP5   303	 /* type of dyad synapse to connect with */
#define AUTAPSE5   304	 /* synapse back to presynaptic node */
#define SYNNUM5    305	 /* number of synapses per connection */
#define SENSCA5    306	 /* synaptic release sensitivity calcium */
#define SRRPOOL5   307	 /* synaptic readily releasable pool */
#define SRRPOOLG5  308	 /* synaptic readily releasable pool gain */
#define SMRRPOOL5  309	 /* synaptic max readily releasable pool */
#define SMAXRATE5  310	 /* maximum sustained synaptic release rate */
#define SGAIN5     311	 /* synaptic gain */
#define SVGAIN5    312	 /* synaptic vgain */
#define SDURH5     313	 /* synaptic high pass time const. */
#define SNFILTH5   314	 /* synaptic high pass nfilt */
#define SHGAIN5    315	 /* synaptic high pass gain */
#define SHOFFS5    316	 /* synaptic high pass offset */
#define SVSIZ5     317	 /* synaptic vesicle size */
#define SCOND5     318	 /* synaptic conductance */
#define SCMUL5     319	 /* synaptic conductance mult for region */
#define SCGRAD5    320	 /* synaptic conductance gradient from soma */
#define SEGRAD5    321	 /* synaptic conductance expon grad fr soma */
#define STHRESH5   322	 /* synaptic threshold */
#define SVNOISE5   323	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV5      324	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR5      325	 /* synaptic event time const. */
#define SFALL5     326	 /* synaptic event fall time const. */
#define SNFILT5    327	 /* synaptic vesicle nfilt */
#define STRCONC5   328	 /* synaptic transmitter concentration. */
#define SRESP5     329	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA5      330	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX5   331	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM5     332	 /* synaptic postsyn Ca pump Km. */
#define SCNFILT5   333	 /* second mesng. nfilt */
#define SCDUR5     334	 /* second mesng. time const. */
#define SCGAIN5    335	 /* synaptic second messenger gain */
#define SCOFF5     336	 /* synaptic second messenger offset */
#define SCNOISE5   337	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN5    338	 /* number of channels */
#define SUNIT5     339	 /* synaptic channel unitary conductace */
#define SVREV5     340	 /* synaptic reversal potential */

#define CELPOST6   341	 /* cell type to connect to (neg, no conn) */
#define CONPOST6   342	 /* connection number for postsyn cell */
#define CELDIV6    343	 /* number of postsyn cells to connect to */
#define GROWPRE6   344	 /* grow when making conn to postsyn cell */
#define SYNREG6    345	 /* synaptic region in presyn dendritic tree */
#define SYNREGP6   346	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC6   347	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI6   348	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO6   349	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI6   350	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO6   351	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG6    352	 /* angle for postsynaptic cell */
#define SYNRNG6    353	 /* range of angles for postsynaptic cell */
#define USEDYAD6   354	 /* synapse is dyad using preexisting type */
#define DYADTYP6   355	 /* type of dyad synapse to connect with */
#define AUTAPSE6   356	 /* synapse back to presynaptic node */
#define SYNNUM6    357	 /* number of synapses per connection */
#define SENSCA6    358	 /* synaptic release sensitivity calcium */
#define SRRPOOL6   359	 /* synaptic readily releasable pool */
#define SRRPOOLG6  360	 /* synaptic readily releasable pool gain */
#define SMRRPOOL6  361	 /* synaptic max readily releasable pool */
#define SMAXRATE6  362	 /* maximum sustained synaptic release rate */
#define SGAIN6     363	 /* synaptic gain */
#define SVGAIN6    364	 /* synaptic vgain */
#define SDURH6     365	 /* synaptic high pass time const. */
#define SNFILTH6   366	 /* synaptic high pass nfilt */
#define SHGAIN6    367	 /* synaptic high pass gain */
#define SHOFFS6    368	 /* synaptic high pass offset */
#define SVSIZ6     369	 /* synaptic vesicle size */
#define SCOND6     370	 /* synaptic conductance */
#define SCMUL6     371	 /* synaptic conductance mult for region */
#define SCGRAD6    372	 /* synaptic conductance gradient from soma */
#define SEGRAD6    373	 /* synaptic conductance expon grad fr soma */
#define STHRESH6   374	 /* synaptic threshold */
#define SVNOISE6   375	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV6      376	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR6      377	 /* synaptic event time const. */
#define SFALL6     378	 /* synaptic event fall time const. */
#define SNFILT6    379	 /* synaptic vesicle nfilt */
#define STRCONC6   380	 /* synaptic transmitter concentration. */
#define SRESP6     381	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA6      382	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX6   383	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM6     384	 /* synaptic postsyn Ca pump Km. */
#define SCNFILT6   385	 /* second mesng. nfilt */
#define SCDUR6     386	 /* second mesng. time const. */
#define SCGAIN6    387	 /* synaptic second messenger gain */
#define SCOFF6     388	 /* synaptic second messenger offset */
#define SCNOISE6   389	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN6    390	 /* number of channels */
#define SUNIT6     391	 /* synaptic channel unitary conductace */
#define SVREV6     392	 /* synaptic reversal potential */

#define CELPOST7   393	 /* cell type to connect to (neg, no conn) */
#define CONPOST7   394	 /* connection number for postsyn cell */
#define CELDIV7    395	 /* number of postsyn cells to connect to */
#define GROWPRE7   396	 /* grow when making conn to postsyn cell */
#define SYNREG7    397	 /* synaptic region in presyn dendritic tree */
#define SYNREGP7   398	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC7   399	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI7   400	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO7   401	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI7   402	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO7   403	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG7    404	 /* angle for postsynaptic cell */
#define SYNRNG7    405	 /* range of angles for postsynaptic cell */
#define USEDYAD7   406	 /* synapse is dyad using preexisting type */
#define DYADTYP7   407	 /* type of dyad synapse to connect with */
#define AUTAPSE7   408	 /* synapse back to presynaptic node */
#define SYNNUM7    409	 /* number of synapses per connection */
#define SENSCA7    410	 /* synaptic release sensitivity calcium */
#define SRRPOOL7   411	 /* synaptic readily releasable pool */
#define SRRPOOLG7  412	 /* synaptic readily releasable pool gain */
#define SMRRPOOL7  413	 /* synaptic max readily releasable pool */
#define SMAXRATE7  414	 /* maximum sustained synaptic release rate */
#define SGAIN7     415	 /* synaptic gain */
#define SVGAIN7    416	 /* synaptic vgain */
#define SDURH7     417	 /* synaptic high pass time const. */
#define SNFILTH7   418	 /* synaptic high pass nfilt */
#define SHGAIN7    419	 /* synaptic high pass gain */
#define SHOFFS7    420	 /* synaptic high pass offset */
#define SVSIZ7     421	 /* synaptic vesicle size */
#define SCOND7     422	 /* synaptic conductance */
#define SCMUL7     423	 /* synaptic conductance mult for region */
#define SCGRAD7    424	 /* synaptic conductance gradient from soma */
#define SEGRAD7    425	 /* synaptic conductance expon grad fr soma */
#define STHRESH7   426	 /* synaptic threshold */
#define SVNOISE7   427	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV7      428	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR7      429	 /* synaptic event time const. */
#define SFALL7     430	 /* synaptic event fall time const. */
#define SNFILT7    431	 /* synaptic vesicle nfilt */
#define STRCONC7   432	 /* synaptic transmitter concentration. */
#define SRESP7     433	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA7      434	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX7   435	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM7     436	 /* synaptic postsyn Ca pump Km. */
#define SCNFILT7   437	 /* second mesng. nfilt */
#define SCDUR7     438	 /* second mesng. time const. */
#define SCGAIN7    439	 /* synaptic second messenger gain */
#define SCOFF7     440	 /* synaptic second messenger offset */
#define SCNOISE7   441	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN7    442	 /* number of channels */
#define SUNIT7     443	 /* synaptic channel unitary conductace */
#define SVREV7     444	 /* synaptic reversal potential */

#define CELPOST8   445	 /* cell type to connect to (neg, no conn) */
#define CONPOST8   446	 /* connection number for postsyn cell */
#define CELDIV8    447	 /* number of postsyn cells to connect to */
#define GROWPRE8   448	 /* grow when making conn to postsyn cell */
#define SYNREG8    449	 /* synaptic region in presyn dendritic tree */
#define SYNREGP8   450	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC8   451	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI8   452	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO8   453	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI8   454	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO8   455	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG8    456	 /* angle for postsynaptic cell */
#define SYNRNG8    457	 /* range of angles for postsynaptic cell */
#define USEDYAD8   458	 /* synapse is dyad using preexisting type */
#define DYADTYP8   459	 /* type of dyad synapse to connect with */
#define AUTAPSE8   460	 /* synapse back to presynaptic node */
#define SYNNUM8    461	 /* number of synapses per connection */
#define SENSCA8    462	 /* synaptic release sensitivity calcium */
#define SRRPOOL8   463	 /* synaptic readily releasable pool */
#define SRRPOOLG8  464	 /* synaptic readily releasable pool gain */
#define SMRRPOOL8  465	 /* synaptic max readily releasable pool */
#define SMAXRATE8  466	 /* maximum sustained synaptic release rate */
#define SGAIN8     467	 /* synaptic gain */
#define SVGAIN8    468	 /* synaptic vgain */
#define SDURH8     469	 /* synaptic high pass time const. */
#define SNFILTH8   470	 /* synaptic high pass nfilt */
#define SHGAIN8    471	 /* synaptic high pass gain */
#define SHOFFS8    472	 /* synaptic high pass offset */
#define SVSIZ8     473	 /* synaptic vesicle size */
#define SCOND8     474	 /* synaptic conductance */
#define SCMUL8     475	 /* synaptic conductance mult for region */
#define SCGRAD8    476	 /* synaptic conductance gradient from soma */
#define SEGRAD8    477	 /* synaptic conductance expon grad fr soma */
#define STHRESH8   478	 /* synaptic threshold */
#define SVNOISE8   479	 /* 1 -> vesicle noise, override, vnoise=0 */
#define SCOV8      480	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR8      481	 /* synaptic event time const. */
#define SFALL8     482	 /* synaptic event fall time const. */
#define SNFILT8    483	 /* synaptic vesicle nfilt */
#define STRCONC8   484	 /* synaptic transmitter concentration. */
#define SRESP8     485	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA8      486	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX8   487	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM8     488	 /* synaptic postsyn Ca pump Km. */
#define SCNFILT8   489	 /* second mesng. nfilt */
#define SCDUR8     490	 /* second mesng. time const. */
#define SCGAIN8    491	 /* synaptic second messenger gain */
#define SCOFF8     492	 /* synaptic second messenger offset */
#define SCNOISE8   493	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN8    494	 /* number of channels */
#define SUNIT8     495	 /* synaptic channel unitary conductace */
#define SVREV8     496	 /* synaptic reversal potential */

#define CELPOST9   497	 /* cell type to connect to (neg, no conn) */
#define CONPOST9   498	 /* connection number for postsyn cell */
#define CELDIV9    499	 /* number of postsyn cells to connect to */
#define GROWPRE9   500	 /* grow when making conn to postsyn cell */
#define SYNREG9    501	 /* synaptic region in presyn dendritic tree */
#define SYNREGP9   502	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC9   503	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI9   504	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO9   505	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI9   506	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO9   507	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG9    508	 /* angle for postsynaptic cell */
#define SYNRNG9    509	 /* range of angles for postsynaptic cell */
#define USEDYAD9   510	 /* synapse is dyad using preexisting type */
#define DYADTYP9   511	 /* type of dyad synapse to connect with */
#define AUTAPSE9   512	 /* synapse back to presynaptic node */
#define SYNNUM9    513	 /* number of synapses per connection */
#define SENSCA9    514	 /* synaptic release sensitivity calcium */
#define SRRPOOL9   515	 /* synaptic readily releasable pool */
#define SRRPOOLG9  516	 /* synaptic readily releasable pool gain */
#define SMRRPOOL9  517	 /* synaptic max readily releasable pool */
#define SMAXRATE9  518	 /* maximum sustained synaptic release rate */
#define SGAIN9     519	 /* synaptic gain */
#define SVGAIN9    520	 /* synaptic vgain */
#define SDURH9     521	 /* synaptic high pass time const. */
#define SNFILTH9   522	 /* synaptic high pass nfilt */
#define SHGAIN9    523	 /* synaptic high pass gain */
#define SHOFFS9    524	 /* synaptic high pass offset */
#define SVSIZ9     525	 /* synaptic vesicle size */
#define SCOND9     526	 /* synaptic conductance */
#define SCMUL9     527	 /* synaptic conductance mult for region */
#define SCGRAD9    528	 /* synaptic conductance gradient from soma */
#define SEGRAD9    529	 /* synaptic conductance expon grad fr soma */
#define STHRESH9   530	 /* synaptic threshold */
#define SVNOISE9   531	 /* 1 -> vesicle noise, override, vnoise=0 */
#define SCOV9      532	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR9      533	 /* synaptic event time const. */
#define SFALL9     534	 /* synaptic event fall time const. */
#define SNFILT9    535	 /* synaptic vesicle nfilt */
#define STRCONC9   536	 /* synaptic transmitter concentration. */
#define SRESP9     537	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA9      538	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX9   539	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM9     540	 /* synaptic postsyn Ca pump Km. */
#define SCNFILT9   541	 /* second mesng. nfilt */
#define SCDUR9     542	 /* second mesng. time const. */
#define SCGAIN9    543	 /* synaptic second messenger gain */
#define SCOFF9     544	 /* synaptic second messenger offset */
#define SCNOISE9   545	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN9    546	 /* number of channels */
#define SUNIT9     547	 /* synaptic channel unitary conductace */
#define SVREV9     548	 /* synaptic reversal potential */

#define NPARAMS    549	 /* number of neuron parameters */

#define CELPRE       0	/* cell type to connect to (neg, no conn) */
#define CONPRE       1	/* connection number of presyn cell */
#define CELCONV      2	/* number of presyn cells to connect to */
#define GROWPOST     3	/* grow when making conn from presyn cell */
#define NCONNP       4	/* number of connection parameters */

#define CELPOST      0	/* cell type to connect to (neg, no conn) */
#define CONPOST      1	/* connection number for postsyn cell */
#define CELDIV       2	/* number of postsyn cells to connect to */
#define GROWPRE      3	/* grow when making conn to postsyn cell */
#define SYNREG       4	/* synaptic region in presyn dendritic tree */
#define SYNREGP      5	/* synaptic region in postsyn dendritic tree */
#define SYNSPAC      6	/* synaptic spacing in presyn dendritic tree */
#define SYNANNI      7	/* inner dia of annulus in presyn dendr tree */
#define SYNANNO      8	/* outer dia of annulus in presyn dendr tree */
#define SYNANPI      9	/* inner dia of annulus in postsyn dend tree */
#define SYNANPO     10	/* outer dia of annulus in postsyn dend tree */
#define SYNANG      11	/* angle for postsynaptic cell */
#define SYNRNG      12	/* range of angles for postsynaptic cell */
#define USEDYAD     13	/* synapse is dyad using preexisting type */
#define DYADTYP     14	/* type of dyad synapse to connect with */
#define AUTAPSE     15	/* synapse back to presynaptic node */
#define SYNNUM      16	/* number of synapses per connection */
#define SENSCA      17	/* synaptic release sensitivity calcium */
#define SRRPOOL     18	/* synaptic readily releasable pool */
#define SRRPOOLG    19	/* synaptic readily releasable pool gain */
#define SMRRPOOL    20	/* synaptic max readily releasable pool */
#define SMAXRATE    21	/* maximum sustained synaptic release rate */
#define SGAIN       22	/* synaptic gain */
#define SVGAIN      23	/* synaptic vgain */
#define SDURH       24	/* synaptic high pass time const. */
#define SNFILTH     25	/* synaptic high pass nfilt */
#define SHGAIN      26	/* synaptic high pass gain */
#define SHOFFS      27	/* synaptic high pass offset */
#define SVSIZ       28	/* synaptic vesicle size */
#define SCOND       29	/* synaptic conductance */
#define SCMUL       30	/* synaptic conductance mult for region */
#define SCGRAD      31	/* synaptic conductance gradient from soma */
#define SEGRAD      32	/* synaptic conductance expon grad fr soma */
#define STHRESH     33	/* synaptic threshold */
#define SVNOISE     34	/* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV        35	/* 1=Poisson, <1->more regular, gamma dist */
#define SDUR        36	/* synaptic event time const. */
#define SFALL       37	/* synaptic event fall time const. */
#define SNFILT      38	/* synaptic vesicle nfilt */
#define STRCONC     39	/* synaptic transmitter concentration. */
#define SRESP       40	/* synaptic response (ampa,gaba,gj,etc. */
#define SPCA        41	/* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX     42	/* synaptic postsyn Ca pump vmax. */
#define SCAKM       43	/* synaptic postsyn Ca pump Km. */
#define SCNFILT     44	/* second mesng. nfilt */
#define SCDUR       45	/* second mesng. time const. */
#define SCGAIN      46	/* synaptic second messenger gain */
#define SCOFF       47	/* synaptic second messenger offset */
#define SCNOISE     48	/* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN      49	/* number of channels */
#define SUNIT       50	/* synaptic channel unitary conductace */
#define SVREV       51	/* synaptic reversal potential */
#define NSYNP       52	/* number of synaptic parameters */

#define NCONNI      10	/* number of input connection cell types */
#define NCONNO       9	/* number of output connection cell types  */

#define XGLUT        1	/* generic glutamate response */
#define XAMPA        2	/* AMPA (type 1) synaptic response */
#define XAMPA1       3	/* AMPA type 1 synaptic response */
#define XAMPA2       4	/* AMPA type 2 synaptic response */
#define XAMPA3       5	/* AMPA type 3 synaptic response */
#define XAMPA4       6	/* AMPA type 4 synaptic response */
#define XAMPA5       7	/* AMPA type 5 synaptic response */
#define XNMDA        8	/* NMDA type 1 synaptic response */
#define XNMDA2       9	/* NMDA type 2 synaptic response */
#define XKAINATE    10	/* Kainate synaptic response */
#define XMGLUR6     11	/* mGluR6 synaptic response */
#define XGABA       12	/* GABA type 1 synaptic response */
#define XGABA1      13	/* GABA type 1 synaptic response */
#define XGABA2      14	/* GABA type 2 synaptic response */
#define XGABA3      15	/* GABA type 3 synaptic response */
#define XGABA4      16	/* GABA type 4 synaptic response */
#define XGLY        17	/* Glycine synaptic response */
#define XGAPJ       18	/* gap junction synaptic response */
#define XDYAD       19	/* dyad synapse (uses other resp type) */
#define NRESPTYPES  20	/* number of synaptic types */

