# Numerical simulation of mud-flows impacting structures

GRECO Massimo<sup>1</sup> https://orcid.org/0000-0002-3806-9797; e-mail: grecom@unina.it 

DI CRISTO Cristiana<sup>1</sup>* https://orcid.org/0000-0001-6578-4502; e-mail: dicristo@unicas.it 

IERVOLINO Michele<sup>2</sup> https://orcid.org/0000-0001-8344-8334; e-mail: michele.iervolino@unicampania.it 

VACCA Andrea<sup>1</sup> https://orcid.org/0000-0002-7170-2005; e-mail: vacca@unina.it 

*Corresponding author 

1 DICEA, Università di Napoli Federico II, Napoli 80125, Italy 

2 DI, Università della Campania Luigi Vanvitelli, Aversa (CE) 81031, Italy 

Citation: Greco M, Di Cristo C, Iervolino M (2019) Numerical simulation of mud-flows impacting structures. Journal of Mountain Science 16(2). https://doi.org/10.1007/s11629-018-5279-5 

© Science Press, Institute of Mountain Hazards and Environment, CAS and Springer-Verlag GmbH Germany, part of Springer Nature 2019 

Abstract: The study of the interaction of mud-flows with obstacles is important to define inundation zones in urban areas and to design the possible structural countermeasures. The paper numerically investigates the impact of a mud-flow on rigid obstacles to evaluate the force acting on them using two different depth-integrated theoretical models, Single-Phase Model (SPM) and Two-Phase Model (TPM), to compare their performance and limits. In the first one the water-sediment mixture is represented as a homogeneous continuum described by a shearthinning power-law rheology. Alternatively, the twophase model proposed by Di Cristo et al in 2016 is used, which separately accounts for the liquid and solid phases. The considered test cases are represented by a 1D landslide flowing on a steep slope impacting on a rigid wall and a 2D mud dam-break flowing on a horizontal bottom in presence of single and multiple rigid obstacles. In the 1D test case, characterized by a very steep slope, the Two-Phase Model predicts the separation between the two phases with a significant longitudinal variation of the solid concentration. In this case the results indicate appreciable differences between the two models in the estimation of both the wave celerity and the magnitude of the impact, with an overestimation of the peak force when using the Single-Phase Model. In the 2D test-cases, where the liquid and solid phases remain mixed, even if the flow fields predicted by the two models present some differences, the essential features of the interaction with the obstacles, along with the maximum impact force, are comparable. 

Keywords: Mud-Flow; Impact force; Two-phase model; Power-law. 

## Introduction

Gravity-driven debris and mud-flows involving multi-phase mixtures typically occur in mountain areas triggered by heavy rains (Hutter et al. 1996; Iverson and Denlinger 2001) and are characterized by the strong interaction between solid and fluid. These phenomena have a very destructive power, producing loss of human lives, damages and modifying topography (Takahashi 2007). In particular, the expression “mud flow” denotes the motion of a highly-concentrated mixture of water and fine sediments. 

The implementation of appropriate risk mitigation strategies is essential to prevent damages and to reduce losses. For this reason, it is crucial to predict the flow evolution, identifying the paths, the runout distances and the inundated areas (Cui et al. 2011), to evaluate the impact on buildings, piles, vegetations, but also to design structural defense measures. Different solutions may be adopted to obstruct landslides or muddebris flows, to reduce their mobility or to deflect their direction. Some examples of these solutions are barriers, breaking mounds, detention basins, channel diversions, deflecting dykes (Hung et al. 1984; Mizuyam 2008; Jóhannesson et al. 2009). In particular, to reduce the flow mobility and to dissipate flow energy, countermeasures include single or multiple structures, such as rigid or flexible barriers, levees, silt dams and baffles (e.g. Ng et al. 2015; Wang et al. 2017). 

To define the inundation zone in urban areas and to select the most appropriate countermeasures, their location and dimension, the study of the complex interaction between debris or mud-flows with obstacles is extremely important. The objective of the present paper is to investigate the performance of two different theoretical models in reproducing the impact of a mud-flow on rigid obstacles. 

The selection of the most adequate model for reproducing the flow of water-solid mixtures is not a simple task and many options are available. In literature, these flows have been mainly investigated using single-phase description of the flowing medium considered as a homogeneous continuum with a non-Newtonian behavior (e.g. Dent and Lang 1983; Liu and Mei 1989; Coussot 1997; Huang and Garcia 1998). Alternatively watersediment flows have been reproduced with quasisingle phase mixture (e.g. Xia et al. 2018) and twophase (e.g. Iverson 1997; Iverson and Denlinger 2001; Pitman and Le 2005; Greco et al. 2012a; Li et al. 2018a,b) models. 

Different models have been proposed describing debris and mud-flows as a homogeneous continuum medium and employing a non-Newtonian behavior to incorporate the effect of particle interactions (Coussot 1994). The adoption of different formulations depends on the mixture characteristics, such as concentration, and on the type of the solid fraction. Some models account for the presence of a yield stress, such as the Bingham model (e.g. Liu and Mei 1989; Imran et al. 2001; Hewitt and Balmforth, 2013) and the Herschel–Bulkley (e.g. Huang and Garcia 1998; Chanson et al. 2006), while the power-law fluid rheology has been proposed for describing more adequately fluids that do not show any appreciable yield stress at low shear rates (Ng and Mei 1994; Hwang et al. 1994; Perazzo and Gratton 2004). Many studies have been performed for investigating the applicability and the characteristics of the power-law fluids (e.g. Burger et al. 2010; Di Cristo et al. 2014; Turnbull et al. 2015; Campomaggiore et al. 2016; Di Cristo et al. 2018a), suggesting that this rheology is suitable for reproducing the behavior of magmas, mining residuals and fine sediment-water mixtures, which are encountered in flows with a finite fraction of sediments such as mud-flows (e.g. Ng and Mei 1994; Sonder et al. 2006; Longo et al. 2015). For instance, a shear-thinning power-law model is the most appropriate in characterizing the rheology of the soil collected in Cervinara (Avellino, South Italy), where a catastrophic landslide occurred in 1999 (Carotenuto et al. 2015) or of the natural estuarine mud dredged from Haihe River in Tianjin and Mazhou Island near Shenzhen (Zhang et al. 2010). 

Iverson (1997) suggested that these singlephase models cannot capture the interactions between the fluid and solid phases that are crucial for the description of the observed behavior of debris flows, recommending the use of two-phase models. However, some of them are essentially quasi single-phase because they neglect the difference between sediment and water velocities; moreover, the role of the pore fluid is parametrically incorporated (Iverson and Denlinger 2001; Pudasaini et al. 2005; Fernandez-Nieto et al. 2008). Without considering the different velocities of the two phases, the drag force is not accounted for. The formulation proposed by Pitman and Le (2005) or its variant suggested by Pelati et al. (2008) accounts for the mass and momentum equations for both the solid and fluid components, including the drag forces, but neglects the viscous effect on fluid phase. Meng and Wang (2016) improved the model of Pitman and Le (2005) considering the contribution of fluid viscosity. 

A significant improvement in two-phase modelling is represented by the model by Pudasaini (2012), which considered three innovative aspects (enhanced non-Newtonian viscous stress, virtual mass and generalized drag forces) for better reproducing important physical characteristics of debris flow dynamic. The model is very complex and not easy to apply, but it showed very good performances in different test cases. Successively, He et al. (2014) obtained results comparable with the ones of Pudasaini (2012) with a model in which the solid stress is modelled through the Mohr-Coulomb plasticity and a Newtonian viscous stress is considered for the liquid. 

Alternatively, Greco et al. (2012a), starting from Di Cristo et al. (2006), proposed a two-phase model, which separately considers the liquid and the solid phases accounting for the difference between their velocities. The model has been applied for reproducing fast transient flows, involving sediment transport in the form of bedload and bottom deformation. Later on, in Di Cristo et al. (2016) and in Di Cristo et al. (2018b), the model has been extended to include also the suspended load and the occurrence of bottom mass failures due to slope instability, respectively. Recently, Li et al. (2018a) presented a two-phase model for reproducing a debris flow over a fixed bed, which incorporates interphase and particleparticle interaction, fluid and solid fluctuations and multiple grain sizes, showing the crucial role of the stresses due to fluctuations and of the adequate estimation of the bed shear stress for reproducing the phenomenon. The same authors extended the model applicability to erodible beds, incorporating the mass exchange between bed and flow with the introduction of a new relationship and adopting a new closure model for estimating the bed shear stress, showing a good agreement with experimental laboratory data Li et al. (2018b). Xia et al. (2018) proposed a quasi-single-phase mixture model in which the stresses due to fluctuations are incorporated and compared its performance with a traditional single-phase mixture and the models by Li et al. (2018a). The single-phase model performs better than the traditional one and even if it the results of the two-phase models are relatively better, it is still competitive in terms of computational cost. 

In structural engineering applications the dynamic impact of a debris flow on a structure is often evaluated by empirical formulas adopting the hydrostatic pressure or the impact velocity of the incident flow multiplied by a safety factor larger than one. However, this practice may produce an uncorrect estimation of the force (Cui et al. 2015; Sovilla et al. 2016). Therefore, to give accurate indications to engineer and structural designers, it is crucial to reproduce the principal kinetic characteristics of the dynamic impact. 

Many literature studies investigated the impact of dry granular flow and snow avalanches on structures (e.g. Tai et al. 2001; Faug 2005; Chiou et al. 2005; Teufelsbauer et al. 2009; Cui and Gray 2013) with special attention to the configuration with multiples obstacles, such as arrays of baffles (e.g., Ng et al.2015). The interaction of viscous debris flows with obstacles and the evaluation of the impact force has been object of both experimental and numerical researches (e.g. Canelli et al. 2012; Scheidl et al. 2013). Cui et al. (2015) carried out laboratory experiments to study the impact force of a viscous debris flow, measuring separately the dynamic pressure of the slurry and the impact pressure of the coarse grains, through a wavelet analysis of the signal. An empirical model is proposed for predicting the slurry impact, while the grain pressure is random and the particles impact frequency in the flow front is larger than in the body. Vagnon and Segalini (2016) performed a set of small scale experiments of a debris flow impacting on a rigid barrier, proposing a new equation for estimating the impact force, which considers flow characteristics, material properties and barrier dimensions. However, owing to the scale effects, the direct application to real cases of results deduced in laboratory may be not straightforward (Iverson 1997; Vagnon and Segalini 2016). 

Among the numerical studies, particular interest plays the study of debris flow impact in presence of multiple obstacles. Recently Gao et al. (2017) simulated the impact pressure of a debris propagating in urban areas with a depth-integrated continuum model, able to consider building blockage effects, bed erosion and deposition changing the mixture concentration. The obstacles increase the depth and velocity of the flow as the debris tends to run up and to deposit in front of the buildings increasing the impact pressure, in line with field observations and experiments. Kattel et al. (2018) modelled a two-phase debris flow as a mixture of a solid particles and viscous fluid down an inclined surface with tetrahedral obstacles of different dimensions, number and orientation, using the quasi-three-dimensional model of Pudasaini (2012). It has been shown that the presence of obstacles may increase the solid and fluid phases separation and it strongly influences the flow spreading, the run out and the deposition. 

Differently of the debris flows, a limited number of studies have been carried out with reference to the interaction of mud-flows with obstacles. Tiberghien et al. (2007) and Laigle and Labbe (2017) studied the impact of mudflows experimentally and numerically, respectively. In the former paper, laboratory experiments have been performed for measuring mud-flow velocity and pressure close to a rigid barrier, demonstrating the existence of two distinct impact regimes associated with supercritical and subcritical incident mudflows. In Laigle and Labbe (2017) the impact of dam-break mud-flow is numerically simulated using the SPH method and the Herschel-Bulkley rheology has been considered. The model, validated on benchmarks, is shown to be able to reproduce the local characteristics of the flow near the obstacle, the length of the dead-zone of fluid at rest which forms upstream of the obstacle and the pressure. 

The present research numerically investigates the interaction of mud-flows on rigid obstacles using depth-integrated models. The depthintegrated schematization, valid when the length scale normal to the bottom is very small compared to longitudinal and transverse length scales, has been widely assumed for simulating the behaviour of different kinds of earth-surface flows, such as dam breaks (e.g., Wu and Wang 2007; Soarez-Frasao et al. 2012), debris/mud events (e.g. O'Brien et al. 1993; Hübl and Steinwendtner 2001; Pudasaini 2012; Iverson and George 2014) and to evaluate associated forces on rigid structures (Shige-eda and J.Akiyama 2003; Bukreev 2009; Kattel et al. 2018). 

In the study two different depth-integrated models have been used. 

The first one essentially consists in the Single-Phase Model (SPM) with the slurry rheology described by a homogeneous power-law fluid. In particular, the model proposed by Ng and Mei (1994), deduced in laminar conditions under the boundary-layer approximation, has been adopted. 

Such a model, although less rigorous than the ones deduced through the asymptotic expansions of solutions of the Cauchy Momentum equations (e.g. Fernandez-Nieto et al. 2010; Noble and Vila 2010), is widely used in environmental applications owing to its simplicity. The validity of laminar flow condition essentially derives from the observation that turbulent non-Newtonian flows occur in rare situations (Rudman et al. 2004; Rudman and Blackburn 2006). Moreover, a physically consistent turbulent model for generalized non-Newtonian fluid is still missing (Gori and Boghi 2011, 2012), despite some progress have been made in the Direct Numerical Simulation (DNS) investigation (Gavrilov and Rudyak 2016, 2107). 

The second one is the two-phase shallow-water model proposed by Di Cristo et al. (2016). In this case the term Two-Phase Model (TPM) is be intended in the sense of Euler-Euler model (Sharma et al. 2017; Wang et al. 2010) and not as Euler-Lagrange model (Kolesnichenko and Shiriaev 2002; Morabito et al. 2004). 

The objective of the present research is to compare performance and limits of the two approaches in reproducing the complex interaction between the mud-flows and the rigid obstacles. Due to the difficulty in deciding which approach is preferable, the main goal of the comparison is to understand if there are significant differences in the simulation results and, in particular, in the impact force estimation. Both one-dimensional and two-dimensional test-cases have been considered and discussed. 

The article is organized as follows. In section 1 both the Single and the Two Phase Models are briefly described. Some details concerning their numerical solution are also given. In Section 2 the application of both models for reproducing the impact of mud-waves, originating by a dam break, with rigid obstacle, is presented. Both 1D and 2D tests have been considered. As far as the latter is concerned, two different configurations, obtained varying the number of rigid obstacles, have been simulated. 

## 1 Model Description

## 1.1 Single-Phase Model (SPM)

The depth-integrated equations of a onedimensional unsteady, gradually-varied, laminar flow of a layer of power-law fluid over a not erodible bed proposed by Ng and Mei (1994), have been straightforwardly extended to the twodimensional case. Denoting with x and y the directions in the horizontal plane and with t the time, the dimensional two-dimensional governing equations (in conservative variables) read: 

$$
\frac {\partial \mathbf {U} _ {S P M}}{\partial t} + \frac {\partial \hat {\mathbf {F}} _ {S P M}}{\partial x} + \frac {\partial \hat {\mathbf {G}} _ {S P M}}{\partial y} = \hat {\mathbf {S}} _ {S P M}\tag{1}
$$

where $\mathbf { U } _ { S P M } = \left[ h _ { m } \quad h _ { m } U _ { m } \quad h _ { m } V _ { m } \right] ^ { T }$ denotes the unknowns vector, in which $h _ { \mathrm { m } }$ is the flow depth and ${ \bf U } _ { m } = \left[ U _ { m } \quad V _ { m } \right] ^ { T }$ is fluid velocity. The expressions of the flux functions $\hat { \mathbf { F } } _ { S P M }$ and $\hat { \mathbf { G } } _ { S P M }$ of the source term $\hat { \mathbf { s } } _ { s P M }$ are: 

$$
\hat {\mathbf {F}} _ {S P M} = \left[ \begin{array}{c} h _ {m} U _ {m} \\ \beta h _ {m} U _ {m} ^ {2} + g \cos \theta_ {x} \frac {h _ {m} ^ {2}}{2} \\ \beta h _ {m} U _ {m} V _ {m} \end{array} \right]
$$

$$
\hat {\mathbf {G}} _ {S P M} = \left[ \begin{array}{c} h _ {m} V _ {m} \\ \beta h _ {m} U _ {m} V _ {m} \\ \beta h _ {m} V _ {m} ^ {2} + g \cos \theta_ {y} \frac {h _ {m} ^ {2}}{2} \end{array} \right]
$$

$$
\hat {\mathbf {S}} _ {S P M} = \left[ \begin{array}{c c} \mathbf {o} & \\ g h _ {m} \sin \theta_ {x} - \frac {\tau_ {m , x}}{\rho_ {m}} \\ g h _ {m} \sin \theta_ {y} - \frac {\tau_ {m , y}}{\rho_ {m}} \end{array} \right]\tag{2}
$$

in which $g$ denotes the gravity acceleration, $\theta _ { x }$ (resp. $\theta _ { y } )$ is the bed slope respect to the horizontal plane in the x (resp. y) direction, considered constant and positive. $\beta = 2 { \left( 2 n + 1 \right) } / { \left( 3 n + 2 \right) }$ is the dimensionless momentum flux correction factor (Di Cristo et al. 2013) and $\rho _ { \mathrm { m } }$ is the density of the mixture, which is assumed constant (Ng and Mei 1994). ${ \pmb \tau _ { \mathrm { m } } }$ represents the bottom shear stress, whose expression, for a power-law fluid, is: 

$$
\pmb {\tau} _ {m} = \mu_ {m} \left(\frac {1 + 2 n}{n} \frac {1}{h _ {m}}\right) ^ {n} \left| \mathbf {U} _ {m} \right| ^ {n - 1} \mathbf {U} _ {m}\tag{3}
$$

It is easy to verify that the system represented by Eq(1) is of hyperbolic type and the expression of the eigenvalues is: 

$$
\begin{array}{c} \lambda_ {1} = \beta \mathbf {U} _ {m} \cdot \mathbf {s} \\ \lambda_ {2, 3} = \beta \mathbf {U} _ {m} \cdot \mathbf {s} \pm \sqrt {\beta (\beta - 1) (\mathbf {U} _ {m} \cdot \mathbf {s}) ^ {2} + g h _ {m}} \end{array} \tag {4}
$$

where s denotes the director cosines of an arbitrary direction in the $( x , y )$ plane. 

System (1) may be solved with any of the numerical schemes commonly employed for Saint-Venant Equations. In the present paper, the Finite Volume solver FIVFLOOD (Leopardi et al. 2002; Greco et al. 2012b) has been adapted. The solver uses a McCormack (predictor-corrector) scheme with a three-point parabolic interpolation of the conserved variables values for evaluating the fluxes at the cells interfaces. The CFL condition, written with reference to the largest eigenvalue (Eq.(4)), has been imposed to define the Δt value. 

Starting from the Single-Phase Model and the corresponding numerical method, the turbulent clear water case has been simulated simply setting the dimensionless momentum flux correction factor equal to one and expressing the bottom shear stress through the Chezy formula. 

## 1.2 Two-Phase Model (TPM)

A reduced version of the two-phase depthintegrated model proposed by Di Cristo et al. (2016) is used here, and it is briefly presented in the following. Differently from the Di Cristo et al. (2016) model, the bed erosion/deposition processes have been neglected and the only bedload dynamics is considered without accounting for the suspended load. Considering a uniformlygraded non-cohesive uniform sediment (with diameter $d )$ and neglecting both lift and virtual (added) mass forces, the dynamics of the mixture is analyzed considering two distinct velocities for sediment and water along with the variability of sediment concentration. The equations are also limited by the absence of inflow and outflow from sidewalls and free-surface. The dimensional governing equations read: 

$$
\frac {\partial \mathbf {U} _ {T P M}}{\partial t} + \frac {\partial \hat {\mathbf {F}} _ {T P M}}{\partial x} + \frac {\partial \hat {\mathbf {G}} _ {T P M}}{\partial y} = \hat {\mathbf {S}} _ {T P M}\tag{5}
$$

in which $\mathbf { U } _ { T P M } = \left[ \delta _ { l } \quad \delta _ { s } \quad \delta _ { l } U _ { l } \quad \delta _ { l } V _ { l } \quad \delta _ { s } U _ { s } \quad \delta _ { s } V _ { s } \right] ^ { T }$ is the conservative unknowns vector, having denoted with $\delta _ { l }$ (resp. $\delta _ { s } )$ the liquid (resp. solid) phase volume for unit bottom surface and with $\mathbf { U } _ { l } = \left[ U _ { l } \quad V _ { l } \right] ^ { T }$ (resp. $\mathbf { U } _ { s } = \left[ U _ { s } \quad V _ { s } \right] ^ { T } \quad )$ the corresponding velocity. 

Denoting with $\rho \boldsymbol { \mathbf { \mathit { l } } }$ (resp. $\rho _ { s } )$ the liquid (resp. solid) density and with r the ratio $r { = } ( \rho _ { s } { - } \rho _ { l } ) / \rho _ { l } ,$ , the flux functions $\hat { \mathbf { F } } _ { T P M }$ and $\hat { \mathbf { G } } _ { T P M }$ and the source term $\hat { \bf s } _ { T P M }$ are given by Di Cristo et al. (2016): 

$$
\hat {\mathbf {F}} _ {T P M} = \left[ \begin{array}{c} \delta_ {l} U _ {l} \\ \delta_ {s} U _ {s} \\ \delta_ {l} U _ {l} ^ {2} + g \cos \theta_ {x} \frac {\left(\delta_ {l} + \delta_ {s}\right) ^ {2}}{2} \\ \delta_ {l} U _ {l} V _ {l} \\ \delta_ {s} U _ {s} ^ {2} + \frac {r}{r + 1} \cos \theta_ {x} \frac {g \delta_ {s} ^ {2}}{2 C _ {s}} \\ \delta_ {s} U _ {s} V _ {s} \end{array} \right]
$$

$$
\hat {\mathbf {G}} _ {T P M} = \left[ \begin{array}{c} \delta_ {l} V _ {l} \\ \delta_ {s} V _ {s} \\ \delta_ {l} U _ {l} V _ {l} \\ \delta_ {l} V _ {l} ^ {2} + g \cos \theta_ {y} \frac {\left(\delta_ {l} + \delta_ {s}\right) ^ {2}}{2} \\ \delta_ {s} U _ {s} V _ {s} \\ \delta_ {s} V _ {s} ^ {2} + \frac {r}{r + 1} \cos \theta_ {y} \frac {g \delta_ {s} ^ {2}}{2 C _ {s}} \end{array} \right]
$$

$$
\hat {\mathbf {S}} _ {T P M} = \left[ \begin{array}{c} 0 \\ 0 \\ g (\delta_ {l} + \delta_ {s}) \sin \theta_ {x} - \left(\frac {\tau_ {l , x}}{\rho_ {l}} + \frac {D _ {x}}{\rho_ {l}}\right) \\ g (\delta_ {l} + \delta_ {s}) \sin \theta_ {y} - \left(\frac {\tau_ {l , y}}{\rho_ {l}} + \frac {D _ {y}}{\rho_ {l}}\right) \\ g \delta_ {s} \frac {r}{r + 1} \sin \theta_ {x} - \left(\frac {\tau_ {s , x}}{\rho_ {s}} - \frac {D _ {x}}{\rho_ {s}}\right) \\ g \delta_ {s} \frac {r}{r + 1} \sin \theta_ {y} - \left(\frac {\tau_ {s , y}}{\rho_ {s}} - \frac {D _ {y}}{\rho_ {s}}\right) \end{array} \right]\tag{6}
$$

In Eqs. (6), D represents the drag force by the water on the solid particles: 

$$
\mathbf {D} = \rho_ {l} C _ {D} \frac {\delta_ {s}}{d} \big (\mathbf {U} _ {l} - \mathbf {U} _ {s} \big) \big | \mathbf {U} _ {l} - \mathbf {U} _ {s} \big |\tag{7}
$$

where $C _ { D }$ is the bulk drag coefficient (Di Cristo et al. 2016). Following Greco et al. (2018), the simplified form of the bottom shear stress acting on the liquid ( <sub>l</sub>) and solid ( <sub>s</sub>) phases are given by 

$$
\pmb {\tau} _ {l} = \rho_ {l} \frac {\mathbf {U} _ {l}}{C _ {\mathrm{Ch}} ^ {2}} \big | \mathbf {U} _ {l} \big | - \pmb {\tau} _ {s}\tag{8}
$$

$$
\pmb {\tau} _ {s} = \mu_ {d} g \delta_ {s} \rho_ {s} \frac {r}{r + 1} \frac {\mathbf {U} _ {s}}{\left| \mathbf {U} _ {s} \right|} + \alpha \rho_ {s} \mathbf {U} _ {s} \left| \mathbf {U} _ {s} \right|\tag{9}
$$

where $\mu _ { \mathrm { d } } ,$ ,  and $C _ { \mathrm { C h } }$ are the dynamic friction, the interparticle collisional stress (Bagnold) and the dimensionless Chezy coefficients, respectively. 

The solid concentration $C _ { \mathrm { s } }$ in Eqs. (6) may be expressed in terms of $\delta _ { \mathrm { s } }$ as it follows: 

$$
C _ {s} = \frac {\delta_ {s}}{K _ {s} d}\tag{10}
$$

where the ratio of the bed-load layer thickness to sediment diameter (K<sub>s</sub>) is given by the following expression (Di Cristo et al. 2016): 

$$
K _ {s} = \frac {k _ {1} ^ {3 / 2}}{\left(1 - k _ {1} \mu_ {s}\right) ^ {3 / 2}} \frac {\left[ \theta_ {c} + \mu_ {s} \delta_ {s} / d \right] ^ {3 / 2}}{\left(\delta_ {s} / d\right) ^ {1 / 2}}\tag{11}
$$

where $\mu _ { \mathrm { s } }$ is the static friction coefficient and $k _ { 1 }$ is a dimensionless coefficient (Di Cristo et al. 2016). 

The eigenvalues of the reduced TPS (Eq.5) are: 

$$
\begin{array}{c} \lambda_ {1} = \mathbf {U} _ {l} \cdot \mathbf {s}; \lambda_ {2} = \mathbf {U} _ {s} \cdot \mathbf {s}; \\ \lambda_ {3, 4} = \mathbf {U} _ {s} \cdot \mathbf {s} \pm \sqrt {\frac {g d r}{2 (r + 1)}} \sqrt {K _ {s} + \frac {\partial K _ {s}}{\partial \delta_ {s}} \delta_ {s}}; \\ \lambda_ {5, 6} = \mathbf {U} _ {l} \cdot \mathbf {s} \pm \sqrt {g (\delta_ {l} + \delta_ {s})} \end{array} \tag {12}
$$

The derivative $\partial K _ { s } / \partial \delta _ { s }$ appearing in Eq. (12) can be easily deduced from Eq. (11). 

The hyperbolic system (5) has been solved using 2D quadrangular meshes, previously developed for analysing 2D conditions even in presence of geofailure (Evangelista et al. 2015; Di Cristo et al. 2018b). The numerical method consists in a mixed cell-centred (CCFV) and node-centred (NCFV) finite-volume. In particular the variables $\delta ,$ $\delta _ { \mathrm { s } }$ U<sub>l</sub> and $\mathbf { U } _ { \mathrm { s } } ,$ are defined at the center of the cell while the bottom elevations, needed to define the bed slope s, are collocated at each node of the cells. The numerical fluxes of Eq. (5) are calculated using the first-order Harten–Lax–Van Leer (HLL) scheme (Harten et al. 1983). 

The time step has been chosen to satisfy the CFL condition, expressed with reference to the largest eigenvalue in Eq. (12). 

## 2 Results

Both single and two-phase models have been applied for analyzing the propagation of dam-break waves over a non-erodible floodplain in presence of rigid obstacles. The test-cases discussed in what follows have been chosen with the aim of representing, among an infinitely wide spectrum of possibilities, two alternative scenarios: one in which slope and friction are expected to dominate the wave propagation, and another one where inertia plays the main role. For the sake of comparison, a third scenario has been simulated in which the flowing medium is the turbulent clearwater, denoted in what follow as Clear Water Model (CWM). 

The impact force exerted against the obstacle is assumed as an integral parameter to concisely highlight the differences resulting from the two alternative approaches, i.e. SPM and TPM. The impact force has been evaluated by numerically computing the following integrals: 

$$
\mathbf {F} _ {S P M} = \rho_ {m} \int_ {\sigma} \left(\beta h _ {m} \mathbf {U} _ {m} U _ {, n} + g \frac {h _ {m} ^ {2}}{2} \hat {\mathbf {n}}\right) \mathrm{d} \sigma\tag{13}
$$

$$
\mathbf {F} _ {T P M} = \rho_ {l} \int_ {\sigma} \left(\delta_ {l} \mathbf {U} _ {l} U _ {l, n} + g \frac {h ^ {2}}{2} \hat {\mathbf {n}}\right) \mathrm{d} \sigma
$$

$$
+ \rho_ {s} \int_ {\sigma} \left(\delta_ {s} \mathbf {U} _ {s} U _ {s, n} + \frac {r}{r + 1} \frac {g \delta_ {s} ^ {2}}{2 C _ {s}} \hat {\mathbf {n}}\right) d \sigma\tag{14}
$$

$$
\mathbf {F} _ {\mathrm{CWM}} = \rho_ {l} \int_ {\sigma} \left(\beta h \mathbf {U} U _ {n} + g \frac {h ^ {2}}{2} \hat {\mathbf {n}}\right) \mathrm{d} \sigma\tag{15}
$$

where $\sigma$ denotes the boundary of the obstacle and nˆ is the corresponding normal unit vector. In Eqs. (13)-(15), the values of the flow variables in the cells adjacent to the obstacles have been used. 

## 2.1 The 1D landslide test-case

The first considered test-case is the collapse of a volume of saturated soil accordingly to the scheme of Figure 1 (Iervolino et al. 2017), inspired to the mudflow occurred in Cervinara (Avellino, South Italy), where a catastrophic landslide occurred in 1999. The fixed volume of slurry $( L _ { m } { = } 1 0 0 \ \mathrm { m } , h _ { 0 } { = } 2 \ \mathrm { m } )$ starts moving downstream in an infinitely wide channel (L=100 m) which ends with a wall. The bed slope $\scriptstyle ( \theta = 2 0 ^ { \circ } )$ is assumed constant. The slurry, represented through the power-law Single-Phase Model with a rheological index $\scriptstyle n = 0 . 0 1$ , is characterized by a solid volumetric concentration of 40% and by a consistency of 71.3 Pa sn. These rheological properties have been described supposing that the fluid in the reservoir consists in a highly-concentrated mixture of water and fine sediment, equivalent to the one characterizing the mud flood occurred in 1999 at Cervinara site (Italy). 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/15f9f7721e960a380ee4922d0c6afdfb3a4d3a0d08b943ccc36798c71ba2e289.jpg)



Figure 1 Sketch of 1D landslide test-case.


In performing the numerical simulation with the Two-Phase Model, the density and the diameter of the sediment particles have been assumed equal to 2.66g/cm3 and 10 m, respectively. The static and dynamic friction angles have been fixed equal to $4 5 ^ { \mathrm { { \circ } } }$ and $3 0 ^ { \circ }$ , respectively; the dimensionless Chezy and the k1 coefficient have been set equal to $C _ { \mathrm { C h } } { = } 2 5$ and 0.66, respectively. In the turbulent clear water scenario, the bottom shear stress is evaluated by Chezy formula with $C _ { \mathrm { C h } } = 2 5$ 

In all the considered situations, the numerical simulations have been performed considering ∆x=0.1 m, and a time-step of ∆t= 1/2048 s which guarantees the fulfillment of the CFL condition for all cases. 

Figure 2a compares the time history of the force (computed with reference to a unitary width in the transversal direction) on the obstacle predicted by the three models. The mudflow dynamics may be discussed by individuating three distinct phases in the time history of the impact force. The first one is the propagation of the landslide until the downstream wall is reached, and it ends when the force suddenly starts to increase. The first impact of the mud wave on the wall represents the second phase, which is characterized by a monotone increase of the wave height and in turn of the impact force up to the attainment of the maximum value. Except for the TPM the flow proceeds as a sequence of waves generated by the upwards reflection from the downstream wall, followed by a further downhill propagation of the slurry against it. The first wave of this sequence is considered as the third phase, which is therefore comprised between the instants when the first and the second maximum force values occur at the downstream wall. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/8ae873467bd26281a9e20ec3cc7889b666c13a70574f7cd77fe9af7550e18499.jpg)



b)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/90b6d78de347fe434abeac77e4df5dc788cbf58c43ba85805a987c8907157a5f.jpg)



t (s)



Figure 2 Time history of computed impact force against the obstacle for (a) Single-Phase Model (SPM), Two-Phase Model (TPM) and Clear Water Model (CWM); (b) Two-Phase Model (TPM), repartition between the liquid and solid phases.


Figure 2a shows that the SPM predicts the shortest time to reach the downstream wall: this is easily explained accounting for the slope-induced acceleration of the slurry, which is far more effective on the highly shear-thinning fluid (powerlaw mixture) compared with both the TPM and the CWM. On the other hand, the TPM exhibits the maximum duration of the approaching phase which exceeds the SPM by about 3 s. Finally, the CWM shows an intermediate behaviour. 

As far as the second phase of the wave dynamics is concerned, the SPM predicts the higher value of the maximum impact force $( F _ { S P M } ^ { \mathrm { m a x } } = 4 . 6 \mathrm { k N } )$ which is nearly twice the value of the TPM $\left( F _ { T P M } ^ { \mathrm { m a x } } = 2 . 4 \mathrm { k N } \right)$ . The TPM value is not very different than the impact force predicted by the CWM $\big ( \mathrm { \it { F } } _ { C W M } ^ { \mathrm { m a x } } = 2 . 7 \mathrm { \ k N \ } \big )$ but slightly smaller. With reference to the two-phase formulation, Figure 2b separately represents the contribution of liquid and solid phases to the resulting force, showing that the latter contributes for only about 20% of the maximum impact force. A close examination of Figure 2b reveals also a small timelag of the impact of the solid phase, which is however almost entirely recovered in the attainment of the maximum force. 

Finally, it is worth of note that the third phase for the TPM does not show the existence of a second maximum. Contrarily to the other two formulations the force does not vary substantially after the first impact. However, a close inspection of Figure 2b reveals a slight but progressive increase of the force due to the sediment phase, which corresponds to the progressive stopover of the sediment up to a stationary rest condition. 

The difference in the maximum values of the impact force between SPM and CWM can be easily explained accounting for both the difference in the medium specific weight ( <sub>CWM</sub>= 9806 $\mathrm { { N / m ^ { 3 } } } ,$ <sub>SPM</sub>= 16317 N/m3) and the higher velocity reached by the SPM wave due to the shear-thinning attitude, which implies higher values of the flow depth at the impact against the wall. 

On the other hand, a more detailed analysis of the wave propagation may be helpful to gain insights on the difference between the impact forces predicted by SPM and TPM. To this aim, Figure 3 shows the computed instantaneous flow depth profiles for the three different formulations at different instants during the first phase of wave propagation. Immediately after the landslide triggering, two waves propagating downstream arise. The first one, originated at x=0 is a “compression” wave, while the second, originated at $\scriptstyle x = - L _ { m } ,$ , is a “rarefaction” one. Owing to the presence of these two waves, the portion of the channel where the flow depth attains the initial value reduces with time. Moreover, the celerity of the rarefaction wave appears larger than that of the compression one. For instance, at t = 6 s Figure 3 indicates that the initial flow depth value is found for $- 3 0 \mathrm { m } < x < 1 0$ m for the TPM (Figure 3c), for -7 $\mathbf { m } < x < 3 3$ m for SPM (Figure 3a) and for -10 m < x $< 3 6$ m for CWM (Figure 3e). 

The SPM behaves similarly to CWM in the rarefaction wave, while the differences are observed on the compression wave, which is far more influenced by the bottom resistance law. The difference in the response of the medium to the slope forcing results in very elongated wave profiles in the downstream region for SPM whereas the downstream tip of the CWM simulation assumes the typical convex shape of the Dressler solution (Dressler 1952). Conversely, the TPM waves are characterised by larger flow depth in the rarefaction region and by a slower and thinner downstream limb, compared with both SPM and CWM counterparts. For instance, at $t { = } 6 \ \mathbf { s } ,$ the flow depth predicted at $x = 5 0$ by SPM m is h<sub>SPM</sub>=1.25 m, close to the clear water value $h _ { C W M } { = } 1 . 4 0 $ , whereas 


a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/36203a5e29f6c03956155d16b59138fee69070a5561c364a5944cfa525102912.jpg)



b)



c)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/03fe3c62735c40b050e1f44f4e58aa63cb649e32a651209f2d0d25c8eaa51382.jpg)



d)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/4825f117c584d969343f694bb33a805927ec6ed225018c6cd2fc060645b5dfee.jpg)



D


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/c3a76ed110d5180244a6648068e4701050135885b2238ba7ba95fbebc4edf8a5.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/250b3cc908028bdf96c8f7559d7cdd3b249d3ee847f60898b58e5b432a8d5b8d.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/6159819ca5d31b60e69ede4548f49d1cd64dc1bf0610947a88f654a4c7ac70be.jpg)



Figure 3 Instantaneous longitudinal flow depth profiles for (a, b) Single-Phase Model (SPM), (c, d) Two-Phase Model (TPM) and (e, f) Clear Water Model (CWM).


b) 

h<sub>TPM</sub>=0.60 m. Furthermore, the compression region of the TPM ends with a bore. The diverse behaviour in the downstream region affects the peak wave heights at the obstacle, depicted in Figure 3: h<sub>TPM</sub>=19.9 m at t=12 s (Figure 3d), whereas at t=10 s h<sub>CWM</sub>=21.8 m (Figure 3f) and $h _ { S P M } { = } 2 0 . 7 \mathrm { m } \left( \mathrm { F i g u r e } 3 \mathrm { b } \right)$ 

The composition of the TPM wave during the approaching and the impact stages can be inferred by Figure 4, showing the repartition of the total flow depth between the two phases. In the downstream part of the wave a progressive separation of the phases is observed, with the sediment wave lagging the liquid one (Figure 4a). This effect is related to the combined action of slope and bottom shear stress, occurring in the considered example. Owing to the different response of the two phases, the actual solid concentration of the mixture at the impact of the wave is smaller than in the initial condition: Figure 4b indicates that at t=12 s only 30% of the total volume is occupied by solid ( <sub>s</sub>=4.7m). 

In conclusion, the observed lower maximum impact force of TPM compared with SPM can be explained by the combination of two factors: a smaller wave height in the TPM and a reduction in the mixture bulk density due to the different dynamics of the two phases. Finally, Figure 4 also shows the occurrence of a further, yet slower, phase separation in the upstream of the wave, where a sediment-only wave progressively develops. 

The above results suggest that in the presence of a significant longitudinal variation of the solid concentration the descriptions provided by the TPM and the SPM may lead to appreciable differences in the estimation of both the wave celerity and the magnitude of the impact force. The latter may be strongly overestimated by using a Single-Phase Model which, owing its nature, cannot correctly account for the change in the concentration of the two phases constituting the mixture along the direction of propagation of the wave. 

## 2.2 The 2D dam-break with obstacles testcase

To study the interaction between mud-flows and single obstacle, in which the inertia may play a main role, the geometric scheme studied in Aureli et al. (2015), is initially considered. The authors experimentally and numerically analyzed the impact of a water wave originated from the sudden removal of the gate on single central rectangular rigid block (Figure 5). The bottom was horizontal and the water level in the reservoir was set equal to 10 cm. In that paper three different mathematical models, i.e. a 2D depth-averaged, a 3D Eulerian and a 3D Smoothed Particle Hydrodynamics (SPH) model, have been considered. A time history of the impact force on the obstacle was measured, by load cells mounted inside the obstacle, and compared with the numerical predictions. It has been shown that once calibrated, all mathematical models are able to reproduce the essential features of the phenomenon. In particular, the 2D shallow water model, although not suitable to accurately reproduce the first few instants of the impact force time history, leads to an error in the peak load estimation within only 10% of the measured values. Therefore, the 2D approach may be considered appropriate for practical applications. Furthermore, the numerical results deduced through the depth integrated model are comparable with the experimental data, as well as with the numerical predictions of far more sophisticated and computationally demanding 3D solvers. 


a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/7b62a374c697909299a8913d9114924978c81de93116a3039028f654b86c5d09.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/db0e8d7ca4adf18bec353fb61e9575b3ef2bb5df9f250097e927d7857480f850.jpg)



Figure 4 Two-Phase Model. Instantaneous longitudinal specific phase volume depth profiles. Solid line: liquid phase; Dashed line: solid phase.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/7a5c86ee654a16aa0ea78ceea496f759790d43ed1940fb31dc51b547bb821dc0.jpg)



Figure 5 Sketch of 2D dam-break with obstacle testcase.


The same geometrical scheme is herein considered with a mud flow (Test 1), reproducing its impact on the block using both the SPM and the TPM. The turbulent clear water case is also reproduced to verify the models’ performances respect the experimental data and the simulation performed by Aureli et al. (2015). The geometry has been successively modified adding two half size obstacles symmetrically on the sides of the central block (Test 2). 

The assumed mud fluid is the same highlyconcentrated mixture of water and fine sediment considered in the previous test-case. Numerical simulations have been performed considering $\Delta x { = } \Delta y { = } 0 . 0 0 5 \mathrm { ~ m } .$ , and a time-step t= 1/4096 s which guarantees the fulfillment of the CFL condition. 

Figure 6a compares the time history of the force on the obstacle predicted by both SPM and TPM. The results referring to the turbulent Clear Water Model (CWM) wave resulting from both present calculation and Aureli et al. (2015), are also reported. For the CWM, the present results are in very good agreement with the ones of Aureli et al. (2015), confirming the quality of the present numerical method. 


a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/6636cd0a267ea4e8546acc3d6f70aef669f868aa9c97c7ca2bdf06f84b35786b.jpg)



b)


Compared with the water wave case, the computed force against the obstacle of the mudflow is notably different for both the considered models, with the peak value that approximately doubles up the CWM one: $F _ { C W M } ^ { \mathrm { m a x } } = 6 . 6 \ : \mathrm { N }$ $F _ { S P M } ^ { \operatorname* { m a x } } = 1 1 . 8 \mathrm { ~ N ~ } , F _ { T P M } ^ { \operatorname* { m a x } } = 1 0 . 8 \mathrm { ~ N ~ }$ . This is partially due to the increase of the specific weight of the mixture respect to the clear-water case. In fact, dividing the peak value of the force by the corresponding specific weight these differences are reduced, even if the predictions relative to the mud-flows exceed the clear-water ones of ~80% and ~60% for the SPM and the TPM, respectively. 

The comparison among the curves in Figure 6a suggests that the impact of the turbulent clearwater wave on the obstacle occurs slightly later with respect to both the SPM and TPM, therefore the propagation celerity of the front of the mudwave appears to be slightly larger than the one pertaining to the turbulent clear-water. The force peak value predicted by the TPM is slightly larger (~10%) than the one deduced through the SPM, while the time at which it occurs appears to be almost independent on the model. 

Figure 6b refers only to the TPM, splitting the total force into the solid and liquid contributions. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/50f7227e987b661e2f6aeaa04c831d9532d41a2c9d5e032413f0cfa80bee1926.jpg)



Figure 6 Time history of computed impact force against the obstacle (Test 1) for (a) Single-Phase Model (SPM), Two-Phase Model (TPM) and Clear Water Model (CWM); (b) repartition between the liquid and solid phase of Two-Phase Model (TPM).


c) 


e)



d)


The solid phase front is characterized by a smaller celerity than the liquid one with a contribution to the total force less important than the one of the liquid phase. 

Aiming to provide a deeper analysis of the differences between the two models, the temporal evolution of the flow depth is discussed in the following. Figures 7 to 10 describe the temporal evolution of the flow depth and the instantaneous streamlines at different instants: t = 0.5 s, 1.0 s, 1. 5 s and 2.0 s. In each figure the results of the SPM (a), TPM (b) in terms of total mud flow depth are shown. For the sake of comparison, the corresponding results considering the Turbulent Clear Water case are reported (Figures 7c, 8c, 9c, 10c). As far as the TPM is concerned, Figures 7d, 8d, 9d, 10d (resp. Figures 7e, 8e, 9e, 10e) report the separate contribution of the liquid (resp. solid) phase, respectively. 

At t=0.5 s (see Figure 7), in both SPM and 


a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/30348168bee43013ff28ec758495240a885789cb5d8fc905703a3c9a14e0202a.jpg)



b)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/d5824eaf67f77346653bea55a733ce2a9a50ae5b7fe1fd60ba7ce62f8ad7e47f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/13377b1eb80ff56ae04637080ae6a24c6e54ffa6c06a65798bb26a4ae5f875af.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/564f17bb4e877bd1f79ffd06a3a2302618c5363dbcf97dae5d7a332c1b95aa97.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/78fc6f89cd8f8a7730661212f805e27920c1bca3807dad406d52bad46f0e0b9f.jpg)



Figure 7 Instantaneous flow depth distribution with streamlines superposed for Test 1 at $\mathrm { t } = 0 . 5$ s for (a) Single-Phase Model (SPM), (b) Two-Phase Model (TPM) and (c) Clear Water Model (CWM), (d) Liquid and (e) solid phase from Two-Phase Model (TPM).



a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/f6415cbc8df90725a6b44b312f8ea0fec38088f4aaab7aca0258d465e31175dc.jpg)



b)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/d0fd36a6aed1866972f66bb17fc5b7bcd56f4f71395b458f8a7dbf8bc63da59b.jpg)



d)



e)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/90d1c1d5a30f8704669b9145a7bc06a91fdd4c267f0c5bbc18ad7669a34f1969.jpg)



c)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/7188f2b8e0fc39406c78fe7ea113d7de7942236c6c98bb3aac43b9d88edbfcf0.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/d14d206115a6b5766ffb39a98d2f81b9dc09ddc4a7bc5865512cb1945f0d89ad.jpg)



Figure 8 Instantaneous flow depth distribution with streamlines superposed for Test 1 at t= 1.0 s for (a) Single-Phase Model (SPM), (b) Two-Phase Model (TPM) and (c) Clear Water Model (CWM), (d) Liquid and (e) solid phase from Two-Phase Model (TPM).



e)



b)



a)



c)


a) 


e)


b) 

c) 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/8afdf45fe1d004f3d0c30da2ed2de01904c0b091359fe7d5ed77186d4f971647.jpg)



d)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/6ce4f51952727a599ea0c48c6756b5e8a30d513938d7d7cab6ac271004c2f0f7.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/01f6cb3e61daabafa6bc7a77d1dc37a74b704691a3d6a48430cbec510c40e51d.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/ad0ef31b062d694cce30e675f0c9edbb388ff62e72746984242e2bcb965030d0.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/576dcf176b3f571a22a42ea6a7ae42a998930ac19c98b3a88abb6e2eb6349d69.jpg)



Figure 9 Instantaneous flow depth distribution with streamlines superposed for Test 1 at t= 1.5 s for (a) Single-Phase Model (SPM), (b) Two-Phase Model (TPM) and (c) Clear Water Model (CWM), (d) Liquid and (e) solid phase from Two-Phase Model (TPM).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/70db577b49f18a3afbe8baf5bb9eda5ec7757881caa350e5a90912fb4e187828.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/1c9116f4745b4a9ec2c4de202e0c411f2a9ebf1d889a41c9f24d3e8049ca6323.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/3523d53e9472fc544e1e481a7d774e394ef9822e2d3fe0bb58346cb5b93b67b9.jpg)



d)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/cabd8912b09701e681c48655037abc3a2dc58992fc48f44481d0dd8a533e96f8.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/3949c62dbdd1ded12ea8522e3198e9d8edcbb1fdbd618571b16ce8fee975814b.jpg)



Figure 10 Instantaneous flow depth distribution with streamlines superposed for Test 1 at t= 2.0 s for (a) Single-Phase Model (SPM), (b) Two-Phase Model (TPM) and (c) Clear Water Model (CWM), (d) Liquid and (e) solid phase from Two-Phase Model (TPM).


TPM, the body of the wave has reached the obstacle and it starts to laterally deflect, splitting in two different streams. The comparison between the two mud flow models puts in evidence that SPM predicts the fastest propagation in the streamwise direction (x), in agreement with the force results shown in Figure 6a. Contrarily to the previous testcase, no phase separation occurs in this example, as clearly indicated by Figures 7d and 7e which show that solid phase strictly follows the liquid one. 

Moreover, even the instantaneous streamlines of the solid and liquid phase appear very similar. The results of the TPM shear some similarities with the CWM. 

After the impact (t=1.0 s, Figure 8), the lateral deflection of the wave caused by the obstacle increases, independently of the fluid and of the model. While the interaction of the split stream with the lateral walls is almost completed in the SPM (Figure 8a), it is just starting with the TPM (Figure 8b). Moreover, in the SPM the wave front has arrived at the downstream end of the plotted region, while in the TPM is still far from it. The comparison between Figure 8b and Figure 8c suggests that presence of the solid phase produces only a slightly reduction of the celerity of the wave in the TPM respect the CWM. The absence of any phase separation in the TPM model is again clearly put in evidence by Figures 8d, 8e. 

At t=1.5 s (Figure 9), in both Mixture and Two-Phases models the two streams leave a large dry region behind the obstacle, which appears to be very regular in the SPM. In contrast the CWM predicts that the two currents surround almost completely the obstacle strongly reducing the extension of the dry zone. Such a difference between the mud-models and the turbulent clear water case is still present at t=2.0 s (Figure 10). At t=2.0 s, both models, similarly to the CWM case, predict the presence of recirculation zones close to the lateral walls. Figures 9d-9e and 10d-10e confirm the absence of any phase separation even at t=1.5s and at t=2.0 s, respectively. 

In the Test 2 only the geometry is changed, adding two half-size obstacles. For a complete comparison with the Test 1, Figure 11a reports the time history of the longitudinal force on the central $\left( F _ { 1 } { = } F _ { 1 , x } \right)$ and lateral $( F _ { 2 , x } )$ obstacle predicted by the SPM and TPM, along with the one produced from the CWM. In this case also the transversal component of the force acting on the lateral obstacle is depicted $( F _ { 2 , y } )$ . Figure 11b refers instead only to the TPM, reporting the total force and the part due to the liquid phase. 

As expected, comparing Figures 6a with 11a, the force acting on the central obstacle is essentially unaffected by the presence of the lateral obstacles. 

About the lateral obstacles, both mud-flow models predict that the horizontal and transversal components start acting on it at about $t = 0 . 8 \ s ,$ without any lag respect to the clear water wave. The longitudinal components of the forces predicted by the two mud-flow models behave similarly, while some differences are observed in the transversal components. 

In detail, both models register essentially at the same time (t=1.9 s) the same peak value of the longitudinal components $F _ { 2 , x S P M } ^ { \mathrm { m a x } } = F _ { 2 , x T P M } ^ { \mathrm { m a x } } = 3 . 5 \mathrm { N }$ which exceeds the clear-water value of ~40%. About the transversal component, the SPM produces a higher peak value respect the TPM, both higher respect to the CWM: $F _ { 2 , y S P M } ^ { \mathrm { m a x } } = 4 . 5 \ : \mathrm { N }$ ， $F _ { 2 , y \ T P M } ^ { \mathrm { m a x } } = 3 . 2 \mathrm { N } \ , \ F _ { 2 , y \ C W M } ^ { \mathrm { m a x } } = 2 . 5 \mathrm { N }$ . Similarly to the Test 1, the differences respect to the clear-water case are essentially due to the increased specific weight of the mixture (results not shown). In the SPM the peak value of the transversal component is registered at about t=1.2 s and then a reduction is observed, while the TPM reaches the maximum at t=1.0 s with only a successive slow decrease. 

Figure 11b confirms that in the TPM the contribution of the solid to the maximum total force is about 40% for both the central and the lateral obstacles. 

Figures 12 to 14 report the temporal evolution of the flow depth and the instantaneous streamlines following the dam break at the instants 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/fb41e2c901d94fd78f4369725e33634d1b3ae8e33c4ee2693173b79aacaf42e4.jpg)



t (s)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/be1f51b7026bcdd3dbca18e9c6f94bbe7f24a9c6b221e82edd6a1de65e83c47c.jpg)



t (s)



Figure 11 Time history of computed impact force against the obstacles (Test 2) for (a) Single-Phase Model (SPM), Two-Phase Model (TPM) and Clear Water Model (CWM); (b) repartition between the liquid and solid phase of Two-Phase Model (TPM). Solid lines: $F _ { 1 } ;$ Dashed lines: $F _ { 2 , x ; }$ Dashed-dotted lines: $F _ { 2 , y } .$



a)



d)



e)



d)


1.0 s, 1. 5 s and 2.0 s as in Figures 8 to 10. The plots at t =0.5 s are not reported because they are identical to the ones of Figure 7, since no interaction with the additional obstacles is observed at that time. 

At t=1 s (Figure 12) the flow interacts with all obstacles independently of the fluid and the rheological description adopted for the mixture. The two streams formed from the central obstacle reach the lateral blocks and they are deflected by them. Higher flow depths start to develop at the front of the central obstacle and on the side of the lateral ones. In all cases there are lateral zones close to the dam and a central area in the back of the central obstacle with no fluid. In the SPM and the TPM there is also small dry area in the back of the lateral blocks, not observed for the CWM. Compared with Test 1, the presence of multiple obstacles does not affect the time evolution of the front. At the considered instant, the SPM predicts a 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/bb30b988d27f7845e2cd3676dd0f2de3f8f3066f972db146ea2919fec1a5a409.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/9068b8140c6c67bc8526c0bd5c7d19ff491b98e10c3f6430504a2941a45d04d3.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/babbb62168669fafbd9883396e8380a0025865f4ccc2582821f4c0243a379b06.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/4645e4a7c69f9e0a394e03dea5708d7a4ea2f6ed63a529dc1b901dfab4d60ed9.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/5f45d96814e398f3612289592665139fdf4d90959853e37a46822ab00842abed.jpg)



Figure 12 Instantaneous flow depth distribution with streamlines superposed for Test 2 at t= 1.0 s for (a) Single Phase Model (SPM), (b) Two-Phase Model (TPM) and (c) Clear Water Model (CWM), (d) Liquid and (e) solid phase from Two-Phase Model (TPM).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/b15a35a677a5571044170f061e2a0e058f8e686b17152a2355054bc506e0f3a5.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/1d79e5e582c0de07b82ee61e5788572479eee306c9f790b4aa0f1ee3d87d2b89.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/79056143d1ebf051cca8cc60995a4890ec0203dc942d25e0fea80cb29eb668eb.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/7931044ed68a0c30ee68d56a3bd135c3fae2af29f457475bdba532861605b3ac.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/d04e35938a042ef1e408aee319a83b4c8f7ebdc34075bb02151bae839d6053c3.jpg)



Figure 13 Instantaneous flow depth distribution with streamlines superposed for Test 2 at t= 1.5 s for (a) Single-Phase Model (SPM), (b) Two-Phase Model (TPM) and (c) Clear Water Model (CWM), (d) Liquid and (e) solid phase from Two-Phase Model (TPM).


a) 

c) 


d)



e)


faster streamwise propagation and a wider expansion in the transversal direction respect to the TPM, but the observed depths close to the obstacles are similar in both models explaining the similar impact force observed in Figure 12a at t=1.0 s. Similarly to the Test 1, in the TPM model no phase separation is observed (Figures 11d and 11e) and the presence of the solid phase produces only a slightly reduction of the celerity of the mud wave respect the CWM (Figure 12b and 12c). 

At t=1.5 s (Figure 13) recirculation zones close to the lateral walls are present independently from the models and the fluid. Some differences between the mud and the clear water flows are evident. In the CWM the dry area behind the central obstacle is more reduced respect both the mud flows. The streamlines of the SPM model after the blocks are straight, very similarly to the case with one obstacle (Figure 9). About the TPM an increase of the depth of the solid phase is observed on side of the lateral blocks. 

At t=2.0 s (Figure 14), the presence of lateral obstacles produces a very large recirculation zone in front of the obstacles, in all the examined cases. The principal difference between the Single-Phase and Two-Phases models pertains to the shape of the dry zone past the central obstacle which is larger in SPM than in TPM. Differently the CWM predicts that the flow surrounds almost completely the central obstacle with a very small dry zone past 

it. 

Owing to the geometry of the chosen 2D testcases, the effects of bottom shear stress are very marginal and there is no forcing due to the slope. The above results, taken collectively, suggest that in such a condition the liquid and solid phases remain strongly coupled during the whole event, precluding any phase separation. Although the flow fields predicted by the two models present some differences, the essential features of the forces acting on the obstacles are similarly predicted by the single-phase and the two-phase models. 

## 3 Conclusion

The present work numerically investigates the impact of a mud-flows on rigid obstacles. Two depth-integrated models based on different theoretical approaches have been employed. A Single-Phase Model (SPM) representing the slurry as a power-law fluid and a Two-Phase Model (TPM), which separately considers the liquid and the solid phases, have been considered. The main goal of the comparison is to understand if there are significant differences in the simulation results and specifically in the impact force estimation. 

Both 1D and 2D test cases have been considered. The 1D test case represents the collapse of a volume of saturated soil moving downstream in a deep slope channel which ends with a wall. In this condition, in which the combined action of slope and bottom shear stress is significant, the TPM shows in the downstream part of the wave a progressive separation of the solid and liquid phases. In terms of impact force, the solid phase contributes for only about 20% of the maximum total value. Appreciable differences in the estimation of both the wave celerity and the magnitude of the impact force are observed between the results of the two models. In particular, the peak value of the force of TPM is lower compared with the one of SPM probably for the combined effect of a smaller wave height and a reduction in the mixture bulk density due to the different dynamics of the solid and liquid phase. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/5bf5ce0e33b87a3ed8375384da04fe48c4040d1fcad5aa00fb0858fd1a4af754.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/39d975204b60e51707657f886ea65a5b500d5e73abe469e69e5d20981ba7b8f3.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/80de4293b448a2cb26d13053aaa762455c2aa2bc526960781ef13f876deacbee.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/e137ebc445c14340e11335ac1d708e01f5a2a7bcd5ed3072dfb0ff61fd84e7a1.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/52f884a7-3131-449c-9080-88349170dfe1/968b5bba1693e96d8181acd5e02d247409c18323acfb85f0c30985fa5da09275.jpg)



Figure 14 Instantaneous flow depth distribution with streamlines superposed for Test 2 at t= 2.0 s for (a) Single-Phase Model (SPM), (b) Two-Phase Model (TPM) and (c) Clear Water Model (CWM), (d) Liquid and (e) solid phase from Two-Phase Model (TPM).


The 2D test case reproduces a mud wave originated from the sudden removal of a gate flowing on a horizontal bottom impacting on a single central rectangular rigid block. The geometry has been successively modified adding two half size obstacles symmetrically placed on both sides of the central block (Test 2). About the single obstacle (Test 1), the TPM does not predict phase separation and the estimated value of the peak force is slightly larger (~10%) than the one deduced through the SPM, while the time at which it occurs appears to be almost independent on the model. The presence of the solid phase produces only a slightly reduction of the celerity of the wave in the TPM model with respect to the turbulent clear water case, but in terms of force its contribute on the peak value is still small. 

With multiple obstacles (Test 2), the force 

## References



Aureli F, Dazzi S, Maranzoni A, et al. (2015) Experimental and numerical evaluation of the force due to the impact of a dambreak wave on a structure. Advances in Water Resources 76: 29- 42. https://doi.org/10.1016/j.advwatres.2014.11.009 





Bukreev VI (2009) Force action of discontinuous waves on a vertical wall. Journal Applied Mechanics and Technical Physics. 50(2): 278-283. https://doi.org/10.1007/s10808-009-0037-7 





Burger J, Haldenwang R, Alderman N (2010) Friction factor-Reynolds number relationship for laminar flow of non-Newtonian fluids in open channels of different cross-sectional shapes. Chemical Engineering Science 85: 3549-3556. https://doi.org/10.1016/j.ces.2010.02.040 





Campomaggiore F, Di Cristo C, Iervolino M, et al. (2016) Development of roll waves in power-law fluids with non-uniform initial conditions. Journal of Hydraulic Research 54(3): 289-306. https://doi.org/10.1080/00221686.2016.1140684 





Canelli L, Ferrero AM, Migliazza M, et al. (2012) Debris flow risk mitigation by means of rigid and flexible barriers – experimental tests and impact analysis. Natural Hazard and Earth System 



acting on the central one is essentially unaffected by the presence of the lateral blocks. The longitudinal component of the force predicted by the two mud-flow models behaves similarly, while some differences are observed in the transversal component with the peak value of the SPM larger (~40%) than the TPM. Similarly to the Test 1, in the TPM model no phase separation is observed and the contribution of the solid to the maximum total force is about 40% for both the central and the lateral obstacles. 

The above results, taken collectively, suggest that in the presence of a phase separation the descriptions provided by the TPM and the SPM may lead to appreciable differences in the estimation of both the wave celerity and the magnitude of the impact force, which is overestimated by using a Single-Phase Model. With a different geometry, characterized in particular by a smaller bed slope, the liquid and solid phases remain strongly mixed. In this case, although the flow fields predicted by the two models present some differences, the essential features of the forces acting on the obstacles are similarly predicted by the two representations of the mixture. 

## Acknowledgements

The work described in the present paper was realized in the framework of the project MISALVA, financed by the Italian Minister of the Environment, Land Protection and Sea. CUP H36C18000970005. 



Science 12: 1693-1699. 





https://doi.org/10.5194/nhess-12-1693-2012 





Carotenuto C, Merola MC, Ãlvarez-Romero M, et al. (2015) Rheology of natural slurries involved in a rapid mudflow with different soil organic carbon content. Colloids and Surfaces A 466: 57-65. https://doi.org/10.1016/j.colsurfa.2014.10.037 





Chanson H, Jarny S, Coussot P (2006) Dam Break Wave of Thixotropic Fluid. Journal of Hydraulic Engineering 132 (3): 280–293. 





https://doi.org/10.1061/(asce)0733-9429(2006)132:3(280) 





Chiou MC, Wang Y, Hutter K (2005) Influence of obstacles on rapid granular flows. Acta Mechanica 175: 195-122. https://doi.org/10.1007/s00707-004-0208-9 





Coussot P (1994) Steady, laminar, flow of concentrated mud suspensions in open channel. Journal of Hydraulic Research 32(4): 535-559. https://doi.org/10.1080/00221686.1994.9728354 





Cui P, Gray JMNT (2013) Gravity driven granular free-surace flow around a circular cylinder. Journal of Fluid Mechanics 720: 314- 337. https://doi.org/10.1017/jfm.2013.42 





https://doi.org/10.1016/j.icheatmasstransfer.2012.02.010 https://doi.org/10.1016/j.icheatmasstransfer.2012.02.010 





Cui P, Hu K, Zhuang J, et al. (2011) Prediction of debris-flow area by combining hydrological and inundation simulation methods. Journal of Mountain Science 8: 1-9. https://doi.org/10.1007/s11629-011-2040-8 





Cui P, Chao Z, Lei Y (2015) Experimental analysis on the impact force of viscous debris flow. Earth Surface Process and Landforms 40:1644-1655. https://doi.org/10.1001/esp.3744 





Dent JD, Lang TE (1983) A biviscous modified Bingham model of snow avalanche motion. Annals Glaciology 4: 42-46. https://doi.org/10.3189/S0260305500005218 





Di Cristo C, Iervolino M, Vacca A (2006) Linear stability analysis of a 1-D model with dynamical description of bed load transport. Journal of Hydraulic Research 44: 480-487. https://doi.org/10.1080/00221686.2006.9521699 





Di Cristo C, Iervolino M, Vacca A (2013) Gravity-driven flow of a shear-thinning power-law fluid over a permeable plane. Applied Mathematical Sciences 7(33-36): 1623-1641. https://doi.org/10.12988/ams.2013.13150 





Di Cristo C, Iervolino M, Vacca A (2014) Simplified wave models applicability to shallow mud flows modeled as power-law fluids. Journal of Mountain Sciences 19: 956-965. https://doi.org/10.1007/s11629-014-3065-6 





Di Cristo C, Greco M, Iervolino M, Leopardi A, Vacca A (2016) Twodimensional two-phase depth-integrated model for transients over mobile bed. Journal of Hydraulic Engineering 142(2), 04015043. 





https://dx.doi.org/10.1061/(ASCE)HY.1943-7900.0001024 





Di Cristo C, Iervolino M, Vacca A (2018a) Applicability of kinematic and diffusive models for mud flows: a steady state analysis. Journal of Hydrology 559: 585-595. https://doi org/10 1016/i ihvdrol 2018 02 016 





Di Cristo C, Evangelista S, Iervolino M, et al. (2018b) Dam-break waves over an erodible embankment: experiments and simulations Journal of Hydraulic Research 56(2): 196-210. https://doi.org/10.1080/00221686.2017.1313322 





Dressler RF (1952) Hydraulic resistance effect upon the dam-break functions. Journal of Research of the National Bureau Standards 49(3): 217-225. 





Evangelista S, Greco M, Iervolino M, et al. (2015) A new algorithm for bank-failure mechanisms in 2D morphodynamic models with unstructured grids. International Journal of Sediment Research 30(4): 382-391. https://doi.org/10.1016/j.ijsrc.2014.11.003 





Faug T (2015) Depth-average analytical solution for free-surface granular flow impacting rigid walls down inclines. Physical Review E 92(6). 





Fernandez-Nieto ED, Bouchut F, Bresch D, et al. (2008) A new Savage-Hutter type model for submarine avalanches and generated tsunami. Journal of Computational Physic 227(16): 7720-7754. https://doi.org/10.1016/j.jcp.2008.04.039 





Fernández-Nieto ED, Noble P, Vila JP (2010) Shallow water equations for non-Newtonian fluids. Journal of Non-Newtonian Fluid Mechanics 165(13-14): 712-732. 





Gao L, Zhang LM, Chen HX (2017) Two dimensional simulation of debris flow impact pressure on buildings. Engineering Geology 226: 236-244. https://doi.org/10.1016/j.enggeo.2017.06.012 





Gavrilov AA, Rudyak VY (2016) Reynolds-averaged modeling of turbulent flows of power-law fluids. Journal of Non-Newtonian Fluid Mechanics 227: 45-55. https://doi.org/10.1016/j.jnnfm.2015.11.006 





Gavrilov AA, Rudyak VY (2017) Direct numerical simulation of the turbulent energy balance and the shear stresses in power-law fluid flows in pipes. Fluid Dynamics 52(3): 363-374. https://doi.org/10.1134/S0015462817030048 





Gori F, Boghi A (2011) Two new differential equations of turbulent dissipation rate and apparent viscosity for non-newtonian fluids. International Communications in Heat and Mass Transfer 38(6): 696-703. 





https://doi.org/10.1016/j.icheatmasstransfer.2011.03.003 





Gori F, Boghi A (2012) A three dimensional exact equation for the turbulent dissipation rate of Generalised Newtonian Fluids. International Communications in Heat and Mass Transfer 39(4): 477-485. 





Greco M, Iervolino M, Leopardi A, et al. (2012a) A Two-Phase Model for Fast Geomorphic Shallow Flows. International Journal of 





Sediment Research 27(4): 409-425. 





Greco M, Iervolino M, Vacca A, et al. (2012b) Two-phase modelling of total sediment load in fast geomorphic transients. River Flow 2012, Proc., Int. Conf. on Fluvial Hydraulics, 1, Colegio de Ingenieros Civiles de Costa Rica (CiC): 643-648. 





Greco M, Iervolino M, Vacca A (2018) Analysis of bedform instability with 1-D two-phase morphodynamical models. Advances in Water Resources 120: 50-64. https://doi.org/10.1016/j.advwatres.2017.07.002 





Harten A, Lax PD, van Leer B (1983) On upstream differencing and Godunov-type schemes for hyperbolic conservation laws. SIAM Review 25(1): 35-61. https://doi.org/10.1137/1025002 





He S, Liu W, Ouyang C, et al. (2014) A two-phase model for numerical simulation of debris flow. Natural Hazard and Earth System Science 2: 2151-2183. https://doi.org/10.5194/nhessd-2-2151-2014 





Hewitt DR, Balmforth NJ (2013) Thixotropic gravity currents. Journal of Fluid Mechanics 727: 56-82. https://doi.org/10.1017/jfm.2013.235 





Hubl J, Steinwendtner H (2001) Two-dimensional simulation of two viscous debris flows in Austria. Physics and Chemistry of the Earth-Part C 26(9): 639-644. https://doi.org/10.1016/S1464-1917(01)00061-7 





Huang X, Garcia MH (1998) A Herschel-Bulkley model for mud flow down a slope. Journal of Fluid Mechanics 374: 305-333. https://doi.org/10.1017/S0022112098002845 





Hung O, Morgan GC, Kelerhals R (1984) Quantitative analysis of debris torrent hazard for design of remedial measures. Canadian Geothecnical Journal 21: 663-677. https://doi.org/10.1139/t84-073 





Hutter C, Svendsen B, Rickenmann D (1996) Debris flow modelling: A review, Continuum Mechanics and Thermodynamics 8(1): 1-35. https://doi.org/10.1007/BF01175749 





Hwang CC, Chen JL, Wang JS, et al. (1994) Linear stability of power law liquid film flowing down an inclined plane. Journal of Physics D: Applied Physics 27: 2297-2301. 





Imran J, Harff P, Parker G (2001) A numerical model of submarine debris flows with graphical user interface. Computers & Geosciences 27(6): 717-729. 





Iervolino M, Carotenuto C, Gisonni C, et al. (2017) Impact Forces of a Supercritical Flow of a Shear Thinning Slurry Against an Obstacle. In: Mikoš M, Casagli N, Yin Y, et al. (eds), Advancing Culture of Living with Landslides. WLF 2017. Springer, https://doi.org/10.1007/978-3-319-53485-5_46 





Iverson RM (1997) The physics of debris flows. Review of Geophysics 35(3): 245-296. https://doi.org/10.1029/97RG00426 





Iverson RM, Denlinger RP (2001) Flow of variably fluidized granular masses across three-dimensional terrain: 1. Coulomb mixture theory. Journal of Geophysical Research, Solid Earth 106(B1): 537-552. https://doi.org/10.1029/2000JB900329 





Iverson RM, George DL (2014) A depth-averaged debris-flow model that includes the effect of evolving dilatancy. I Physical basis. Proceeding of Royal Society A Mathematical Physical and Engineering Sciences 470: 1-31. 





Jóhannesson T, Gauer P, Issler P, et al. (2009) The design of avalanche protection dams-Recent practical and theoretical developments. Project Report EUR23339. Climate Change and Natural Hazard Research Area. Series2. European Commission (Available online at: https://hal.archives-ouvertes.fr/hal-00575782/) 





Kattel P, Kafle J, Fischer JT, et al. (2018) Interaction of two-phase debris with obstacles. Engineering Geology 242: 197-217. https://doi.org/10.1016/j.enggeo.2018.05.023 





Kolesnichenko O, Shiriaev AS (2002) Partial stabilization of underactuated Euler–Lagrange systems via a class of feedback transformations. Systems & Control Letters 45(2): 121-132. https://doi.org/10.1016/S0167-6911(01)00170-0 





Laigle D, Labbe M (2017) SPH-based numerical study of the impact of mudflows on obstacles. International Journal of Erosion Control Engineering 10(1): 56-65. https://doi.org/10.1007/s10596-007-9053-y 





https://doi.org/10.1016/j.coldregions.2016.03.002 





Leopardi A, Oliveri E, Greco M (2002) Two-dimensional modeling of flood to map risk prone areas. Journal of Water Resources Planning and Management 128(3): 168-178. https://doi.org/10.1061/(ASCE)0733-9496(2002)128:3(168) 





Li J, Cao ZX, Hu KH, et al. (2018a) A depth-averaged two-phase model for debris flows over fixed beds. International Journal of Sediment Research 33(4): 462-477. https://doi.org/10.1016/j.ijsrc.2017.06.003 





Li J, Cao ZX, Hu KH, et al. (2018b) A depth-averaged two-phase model for debris flows over erodible beds. Earth Surface Processes and Landform 43(4): 817-839. https://doi.org/10.1002/esp.4283 





Liu KF, Mei CC (1989). Slow spreading of a sheet of Bingham fluid on an inclined plane. Journal of Fluid Mechanics 207: 505-529. https://doi.org/10.1017/S0022112089002685 





Longo S, Di Federico V, Chiapponi L (2015) Non-Newtonian powerlaw gravity currents propagating in confining boundaries. Environmental Fluid Mechanics 15: 515. https://doi.org/10.1007/s10652-014-9369-9 





Meng X, Wang Y (2016) Modelling and numerical simulation of twophase debris flows. Acta Geotechnica 11: 1027-1045. https://doi.org/ 10.1007/s11440-015-0418-4 





Mizuyama T (2008) Structural Countermeasures for debris flow disaster. International Journal of Erosion Control Engineering 1(2): 38-43. https://dx.doi.org/10.13101/ijece.1.38 





Morabito F, Teel AR, Zaccarian L (2004) Nonlinear anti-wind-up applied to Euler-Lagrange systems. IEEE Transactions on Robotics and Automation 20(3): 526-537. https://dx.doi.org/10.1109/TRA.2004.824933 





Ng C, Mei CC (1994) Roll waves on a shallow layer of mud modeled as a power-law fluid. Journal of Fluid Mechanics 263: 151-184. https://doi.org/10.1017/S0022112094004064 





Ng CWW, Choi CE, Song D, et al. (2015) Physical modelling of baffled influence on landslide debris mobility. Landslide 12(1): 1- 18. https://doi.org/10.1007/s10346-014-0476-y 





Noble P, Vila JP (2013) Thin power-law film flow down an inclined plane: consistent shallow-water models and stability under largescale perturbations. Journal of Fluid Mechanics 735: 29-60. https://doi.org/10.1017/jfm.2013.454 





O’Brien JS, Julien PY, Fullerton WT (1993) Two-dimensional water flood and mudflow simulation. Journal of Hydraulic Engineering 119(2): 244-261. 





Pelati M, Bouchut F, Mangeney A (2008) A Roe-type scheme for two-phase shallow granular flows over variable topography. ESAIM Mathematical Modelling and Numerical Analysis, 42: 851- 885. https://doi.org/10.1051/m2an:2008029 





Perazzo CA, Gratton J (2004) Steady and traveling flows of a power–law liquid over an incline. Journal of Non-Newtonian Fluid Mechanics 118: 57-64. 





Pitman EB, Le L (2005) A two-fluid model for avalanche and debris flows. Philosophical Transactions of the Royal Society A. Mathematical Physical and Engineering Sciences 363: 1573-1601. https://doi.org/10.1098/rsta.2005.1596 





Pudasaini SP, Wang Y, Hutter K (2005) Modelling debris flows down general channels. Natural Hazard and Earth System Science 5(6): 799-819. https://doi.org/10.5194/nhess-5-799-2005 





Pudasaini SP (2012) A general two-phase debris flow model. Journal of Geophysical Research 117 F03010. https://doi.org/10.1029/2011JF002186 





Rudman M, Blackburn HM, Graham LJW, et al. (2004) Turbulent pipe flow of shear-thinning fluids. Journal of Non-newtonian Fluid Mechanics 118(1): 33-48. 





Rudman M, Blackburn HM (2006) Direct numerical simulation of turbulent non-Newtonian flow using a spectral element method. 





Applied Mathematical Modelling 30(11): 1229-1248. https://doi.org/10.1016/j.apm.2006.03.005 





Sharma R, May J, Alobaid F, et al. (2017) Euler-Euler CFD simulation of the fuel reactor of a 1 MWth chemical-looping pilot plant: Influence of the drag models and specularity coefficient. Fuel 200: 435-446. https://doi.org/10.1016/j.fuel.2017.03.076 





Scheidl C, Chiari M, Kaitna R, et al. (2013) Analysing debris-flow impact model, based on small scale modelling approach. Survey of Geophysics 34(1): 121-140. https://doi.org/10.1007/s10712-012-9199-6 





Shige-eda M, Akiyama J (2003) Numerical and experimental study on two-dimensional flood flows with and without structures. Journal of Hydraulic Engineering 129(10): 817-821. https://doi.org/10.1061/(ASCE)0733-9429(2003)129:10(817) 





Soares-Frazão S, Canelas R, Cao Z, et al. (2012) Dam-break flows over mobile beds: experiments and benchmark tests for numerical models. Journal of Hydraulic Research 50(4): 364-375. https://doi.org/10.1080/00221686.2012.689682 





Sonder I, Zimanowski B, Buttner R (2006) Non-Newtonian viscosity of basaltic magma. Geophysical Research Letter 33: L02303. https://dx.doi.org/10.1029/2005GL024240 





Sovilla B, Faug T, Kohler A, et al. (2016) Gravitational wet avalanche pressure on pylon-like structures. Cold Regions Science and Technology 126: 66-75. 





Tai YC, Gray JMN, Hutter C, et al. (2001) Flow dense avalanches past obstructions. Annals of Glaciology 32: 281-284 https://doi.org/10.3189/172756401781819166 





Takahashi T (2007) Debris Flow: Mechanics, Prediction and Countermeasures. Taylor and Francis, New York, USA. 





Teufelsbauer H, Wang Y, Chou C, et al. (2009) Flow obstacleinteraction in rapid granular avalanches: DEM simulation and comparison with experiments. Granular Matter 11(4): 209-220. https://doi.org/10.1007/s10035-009-0142-6 





Tiberghien D, Laigle D, Naaim M, et al. (2007) Experimental investigation of interaction between mudflow and on obstacle. Proceeding of the International Conference on Debris-Flow Hazard Mitigation: Mechanics, Prediction and Assessment, Chengdu. China. pp 281-292. 





Turnbull B, Bowman ET, McElwaine JN (2015) Debris flows: experiments and modelling. Comptes Rendus Physique 16(1): 86- 96. 





Vagon F, Segalini A (2016) Debris flow impact estimation on a rigid barrier. Natural Hazard and Earth System Science 16: 1691-1697. https://doi.org/10.5194/nhess-16-1691-2016 





Wang Y, Williams KC, Jones MG, et al. (2010) CFD simulation of gas-solid flow in dense phase bypass pneumatic conveying using the Euler-Euler model. Applied Mechanics and Materials 26-28: 1190-1194. 





https://doi.org/10.4028/www.scientific.net/AMM.26-28.1190 





Wang F, Chen X, Chen J, et al. (2017) Experimental study on a debris-flow drainage channel with different types of energy dissipation baffles. Engineering Geology 220: 43-51. https://doi.org/10.1016/j.enggeo.2017.01.014 





Wu W, Wang SS-Y (2007) One dimensional modeling of dam-break flow over movable beds. Journal of Hydraulic Engineering 133(1): 48-58. https://doi.org/10.1061/(ASCE)0733-9429(2007)133:1(48) 





Xia CC, Li J, Cao ZX, et al. (2018) A quasi single-phase model for debris flows and its comparison with a two-phase model. Journal of Mountain Science 15(5): 1071-1089. https://doi.org/10.1007/s11629-018-4886-5 





Zhang X, Bai Y, Ng CO (2010) Rheological Properties of Some Marine Muds Dredged from China Coasts. Proceedings of the 28 International Offshore and Polar Engineering Conference, Beijing, China. pp 455-461. 

