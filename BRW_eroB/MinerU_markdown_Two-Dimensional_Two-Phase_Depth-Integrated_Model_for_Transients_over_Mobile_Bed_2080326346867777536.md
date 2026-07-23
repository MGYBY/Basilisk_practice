# Two-Dimensional Two-Phase Depth-Integrated Model for Transients over Mobile Bed

Cristiana Di Cristo<sup>1</sup>; Massimo Greco<sup>2</sup>; Michele Iervolino<sup>3</sup>; Angelo Leopardi<sup>4</sup>; and Andrea Vacca<sup>5</sup> 

Abstract: Fast geomorphic transients may involve complex scenarios of sediment transport, occurring near the bottom as bed load (i.e., saltating, sliding, and rolling) or as suspended load in the upper portion of the flow. The two sediment transport modalities may even coexist or alternate each other during the same event, especially when the shear stress varies considerably. Modeling these processes is therefore a challenging task, for which the usual representation of the flow as a mixture may result in being unsatisfactory. In the present paper, a new two-phase depth-averaged model is presented that accounts for variable sediment concentration in both bed and suspended loads. Distinct phase velocities are considered for bed load, whereas the slip velocity between the two phases is neglected in the suspended load It is shown that the resulting two-phase model is hyperbolic, and the analytical expression of the eigenvalues is provided. The entrainment/ deposition of sediment between the bottom and the bed load layer is based on a modified van Rijn transport parameter, whereas for the suspended sediment a first-order exchange law is considered. A numerical finite-volume method is used for the simulation of three dam break experiments found in the literature, which are effectively reproduced in terms of both free surface elevation and bottom deformation, confirming the key role played by the solid concentration variability even for two-phase models. DOI: 10.1061/(ASCE)HY.1943-7900 .0001024. This work is made available under the terms of the Creative Commons Attribution 4.0 International license, http:// creativecommons.org/licenses/by/4.0/. 

Author keywords: Two-phase depth-integrated model; Variable concentration; Bed load; Suspended load; Finite-volume method. 

## Introduction

Morphological evolution in river, estuarine, and tidal environments involves the interplay of fluid flow, sediment transport, and loose bed deformation. During extreme events, such as flash floods, avalanche-induced flood waves, debris flows, or dam collapses, the above processes may evolve with comparable time scales. The resulting morphological evolution may lead to dramatic consequences in terms of damages and loss of human lives (Brooks and Lawrence 1999). Analysis and prediction of these fast morphological transients are therefore mandatory for hazard assessment (Sturm 2013). The present paper aims to contribute to this field by presenting a two-phase depth-integrated model suitable for fast unsteady flows, involving sediment transport and bed deformation. 

During unsteady morphological processes, the sediment entrained from the bed is transported through bed load and suspended load. The former occurs under moderate bottom shear stress, whereas the latter pertains to higher bottom shear stress. 

Bed load motion is strongly affected by particle-bottom and particle-particle collisions and by the drag received by the fluid. The suspended load is mainly characterized by the convection by the carrying fluid, often with negligible slip velocity and particle contact. In the presence of a strong spatial and/or temporal variability of the bed shear stress, the two transport modalities may coexist or alternate each other. 

Experimental modeling of fast geomorphic transients encounters strong difficulties. In fact, high-resolution measurements in both time and space of flow field, sediment transport, and bottom deformation are tremendously expensive and beyond the capabilities of most laboratories. With the growing availability of computational resources, the mathematical modeling of these processes is becoming an increasingly interesting alternative for practitioners and researchers. 

The present study follows a deterministic approach, describing the main features of the sediment transport in terms of timeaveraged flow properties only. This approach has the great advantage that the sediment dynamics may be analyzed without a detailed knowledge of the whole process at the price of losing some information concerning the turbulence dynamics. Although this approach is the most used in engineering applications, different analyses have been alternatively developed accounting for the turbulence effect on the sediment transport. For instance, starting from experimental evidence and following a stochastic approach, Papanicolaou et al. (2002a) developed a theoretical model for the inception of sediment motion, accounting for near-bed turbulent structures and bed microtopography. Wu and Chou (2003), incorporating the probabilistic features of the turbulent fluctuations and of the bed-grain geometry, investigated the probability of rolling and lifting for the sediment entrainment. Cheng (2006) showed that the mobility probability of a bed particle may be either enhanced or weakened by an increase of the shear stress fluctuation. In the case of low sediment entrainment, the mobility probability is increased by the turbulence, whereas it is reduced by the shear stress fluctuation if the average bed shear stress becomes relatively high. Wong et al. (2007) designed a detailed experiment to predict the probability density function for the particle virtual velocity and the thickness of the active layer, showing that the statistics of tracer displacements can be related to the macroscopic aspects of the bed load. Furbish et al. (2012) provided a probabilistic definition of the bed load sediment flux. Their formulation is shown to be consistent with experimental measurements and simulations of particle motion. Additionally, either numerical solution of the Reynoldsaveraged Navier-Stokes (e.g., Duran et al. 2012; Marsooli and Wu 2015) equations or of the direct and large eddy simulations (e.g., Keylock et al. 2005; Soldati and Marchioli 2012) of the turbulent flows coupled with sediment particle motion provided useful insights into the role of the coherent structures on erosion/ deposition dynamics. 

In the following only depth-integrated models are considered and discussed. These models do not explicitly account for the dynamics of the very near-bed zone, i.e., the roughness layer. In such a layer, since the flow around sediment particles is strongly three-dimensional and influenced by wakes shed by grains, the velocity profile can significantly deviate from the logarithmic one (Byrd and Furbish 2000; Wohl and Thompson 2000). Considering that the mixing from wakes shed by particles induces a change in the eddy viscosity (Lopez and Garcia 1996; Nikora and Goring 2000; Defina and Bixio 2005). Lamb et al. (2008) assumed a mixing length proportional to the roughness height and derived a parabolic velocity profile, instead of a logarithmic one, in the layer. 

Depth-integrated models may be further distinguished between coupled and decoupled ones. In the coupled models it is assumed that the sediment transport and the bottom evolution develop synchronously (Cao and Carling 2002). On the other hand, decoupled models assume a time-scale hierarchy by which hydrodynamics is usually considered to be faster than the sediment transport and the bottom evolution. 

Common examples of decoupled models are those built up by supplementing a proper fixed-bed hydrodynamic model with a sediment continuity equation (the so-called Exner equation). In the simplest formulation (Graf 1998), the solid discharge is further assumed to instantaneously adapt to the transport capacity, which is estimated by means of empirical relationships proposed for uniform flow conditions (Graf 1998; Wang and Wu 2005). In many real situations this hierarchy is not acceptable, and the application of these models is questionable. Limitations of the decoupled approach have been discussed in the literature (Cao et al. 2002; Garegnani et al. 2011) along with the drawbacks of models based on immediate adaptation of the solid discharge to the transport capacity (Simpson and Castelltort 2006; Di Cristo et al. 2006; Singh et al. 2004; Xia et al. 2010). 

Among the existing coupled (i.e., nonequilibrium) morphological models, a further distinction arises from the representation of the fluid-sediment motion. They may be classified either as mixture or two-phase models, which is the type used herein. To highlight the features of two-phase models, it is useful to first discuss the more popular mixture models. For relatively low solid concentrations, the rheological behavior of the mixture may be represented through clear-water friction law (Wu 2007; Wu and Wang 2007; Sabbagh-Yazdi and Jamshidi 2013). As far as hyperconcentrated mud flows are considered, non-Newtonian constitutive relations able to describe the shear-thinning behavior of the flow are used in the model based on full (Ancey 2012) or simplified (Di Cristo et al. 2014b, c, 2015) wave dynamics. 

The description of a stratified flow with clear water above the mixture leads to the two-layer models, with equal (Fraccarollo and Capart 2002) or distinct (Capart and Young 2002; Li et al. 2013) velocities in the layers. However, within the transport layer no distinction is made between the motion regime of sediments and water. The interaction between mixture and clear-water layers is expressed through an interface shear stress based on the analogy with the multilayer shallow water models. Furthermore, most of these models (Capart and Young 2002; Savary and Zech 2007; Swartenbroekx et al. 2013) assume constant sediment concentration in the transport layer. These models are effective in the analysis of fast morphological transients (Spinewine et al. 2007; Chen and Peng 2006), but the assumption of constant concentration under highly unsteady conditions has been recently questioned. Li et al. (2013) suggested that sediment concentration has to be considered as one of the unknowns of the numerical model, proposing an enhanced two-layer formulation through the application of the fundamental mass conservation law for sediment. Their numerical tests support the conclusion that bed load concentration variability has to be taken into account, if a detailed description of the sediment routing is sought. The mixture models lack any explicit representation of the features of different transport regimes, i.e., bed load and suspended load, which are comprehensively lumped with the behavior of the mixture layer. Furthermore, in these models a hyperbolicity loss may occur in both subcritical and supercritical flow regimes (Savary and Zech 2007; Greco et al. 2008b; Savary and Zech 2008). 

Two-phase modeling is an effective alternative for analyzing the morphohydrodynamics of rivers, debris flows, and snow avalanches (Armanini 2013). Usually, these models are deduced by averaging the conservation principles of mass and momentum for the liquid-solid mixture considered as an equivalent continuous fluid characterized by unique physical characteristics and a unique velocity value, obtaining a phase-averaged system of equations with an unknown variable concentration (e.g., Dewals et al. 2011; Canelas et al. 2013). The system of partial differential equations is hyperbolic, and it may be solved through standard finite volume schemes (Garegnani et al. 2011; Rosatti and Begnudelli 2013). Alternatively, Greco et al. (2012a) proposed a two-phase model that separately considers the liquid and solid phases, accounting for the difference between their velocities and preserving the hyperbolic nature of the system (Evangelista et al. 2013). However, in Greco et al. (2012a) the hypothesis of a constant bed load concentration has been assumed and the suspended load has not been considered. Recent research suggests that these two assumptions should be reconsidered. Indeed, the results by Li et al. (2013), even if referred to mixture models, suggest that the hypothesis of constant bed load concentration may represent a strong limitation. On the other hand, Zhang et al. (2013) recommend that the simulation of both bed load and suspended load may be required to analyze transients with a wide range of shear stress. 

In the present paper a two-phase depth-integrated model is proposed, which is an extension of the preliminary version presented at the River Flow International Conference (Di Cristo et al. 2014a). The model accounts for both the bed and suspended load. As far as the former is concerned, both the liquid-solid velocities difference and the concentration variability are considered. The suspended load is still described assuming the concentration variability, but neglecting the slip velocity between the two phases. The entrainment\deposition of sediments between the bottom and the bed load is evaluated by a formula based on the modified van Rijn mobility parameter, whereas a diffusive vertical flux is assumed to drive the sediments toward the upper region of flow, where the suspended sediment transport occurs. The model is numerically integrated using a finite volume method, and its performance is tested against literature experimental test cases, reporting also the comparison with other existing models. 

The paper is structured in the following way: the proposed model is presented in the next section. In the first subsection the governing equations are given, whereas the closures, the model mathematical characterization, i.e., its hyperbolic nature and a concise presentation of the numerical model are reported in the last two subsections. Then, the results of the model in reproducing experimental data are presented, along with a comparison with other literature models. Finally, the conclusions are drawn. 

## Two-Phase Model

## Governing Equations

In the proposed two-phase model the following hypotheses are assumed: 

• The liquid $( \rho _ { l } )$ and solid $( \rho _ { s } )$ densities are constant; 

• The sediment is uniformly graded (with diameter d) and noncohesive; 

• There is no inflow/outflow from sidewalls and free-surface; and • Standing bed is saturated with a porosity $p .$ 

In the depth-integrated framework, the following shallow-water assumptions are also considered: the vertical components of both acceleration and velocity are neglected; the hydrostatic pressure distribution along the vertical axis is assumed. While these conditions are not strictly verified in the near field of fast geomorphic transients (e.g., during the first instants and in the tip region of a dam break), shallow-water depth-integrated models are widely applied for simulating such events (e.g., Soares-Frazão et al. 2012; Li et al. 2013). In addition, it is supposed that the volume concentration, $C _ { s , b }$ , along the vertical axis of the bed load region is constant and that the suspended sediment passively follows the motion of the fluid phase (Greco et al. 2012b). 

The bed load dynamics is described considering separately the liquid and solid phase, with distinct velocities and accounting for the momentum exchange between them, instead of assuming an equivalent homogeneous fluid with a unique velocity value, i.e., as a water-sediment mixture. Similarly to most of the geophysical flow models (e.g., Pitman and Le 2005; Pudasaini et al. 2005; Pelanti et al. 2008), the lift and virtual (added) mass forces are neglected. As far as the latter force is concerned, Pudasaini (2012) has shown that its introduction in a two-phase model produces a strong coupling (in both time and space) between the streamwise and cross-stream velocity components in the differential terms. However, the inclusion of this force allows only a slight improvement of the model performance in predicting fast processes. On the other side, it has been shown that this additional term, modifying the differential structure of the model, may cause a loss of hyperbolicity and therefore the mathematical well posedness of the system equations is not guaranteed. 

The governing equations, reported in the following, derive from the mass and momentum conservation for the liquid phase [Eqs. (1) and (4)] and solid phase, which moves as bed load [Eqs. (2) and (5)]. Eq. (3) represents the mass conservation of sediment moving as suspended load. Since it is assumed that the sediment velocity is equal to the liquid one in the region where suspended transport occurs, there is no drag between the two phases and therefore the momentum conservation equation for the suspended sediment is not needed. Finally, Eq. (6) is the equation for predicting bed deformation. The complete set of equations reads 

$$
\frac {\partial \delta_ {l}}{\partial t} + \nabla \cdot (\delta_ {l} \mathbf {U} _ {l}) - p e _ {B} = 0\tag{1}
$$

$$
\frac {\partial \delta_ {s , b}}{\partial t} + \nabla \cdot (\delta_ {s, b} \mathbf {U} _ {s}) - (1 - p) e _ {B} + e _ {s, b - s} = 0\tag{2}
$$

$$
\frac {\partial \delta_ {s , s}}{\partial t} + \nabla \cdot (\delta_ {s, s} \mathbf {U} _ {l}) - e _ {s, b - s} = 0\tag{3}
$$

$$
\frac {\partial \delta_ {l} \mathbf {U} _ {l}}{\partial t} + \nabla \cdot (\delta_ {l} \mathbf {U} _ {l} \mathbf {U} _ {l}) + \nabla \left(\frac {g h ^ {2}}{2}\right) + g h \nabla (z _ {B}) + \mathbf {S} _ {l} = 0\tag{4}
$$

$$
\frac {\partial \delta_ {s , b} \mathbf {U} _ {s}}{\partial t} + \nabla \cdot (\delta_ {s, b} \mathbf {U} _ {s} \mathbf {U} _ {s}) + \frac {r}{r + 1} \nabla \left(\frac {g \delta_ {s , b} ^ {2}}{2 C _ {s , b}}\right)
$$

$$
+ g \delta_ {s, b} \frac {r}{r + 1} \nabla (z _ {B}) + \mathbf {S} _ {s b} = 0\tag{5}
$$

$$
\frac {\partial z _ {B}}{\partial t} + e _ {B} = 0\tag{6}
$$

in which t = time; g = gravity acceleration; $r = ( { \rho } _ { s } - { \rho } _ { l } ) / { \rho } _ { l } ;$ and $h = z _ { w } - z _ { B }$ , where $z _ { w }$ and $z _ { B }$ are the free surface and bottom elevation, respectively. In Eqs. $( 1 ) \mathrm { - } ( 5 ) , \delta _ { l }$ denotes the liquid-phase volume for unit bottom surface, $\delta _ { s , b }$ (resp. $\delta _ { s , s } )$ is the solid-phase volume transported as bed (resp. suspended) load for unit bottom surface so that $h = \delta _ { l } + \delta _ { s , b } + \delta _ { s , s } . \mathbf { U } _ { l }$ (resp. $\mathbf { U } _ { s } )$ is the phase-averaged water (resp. solid) velocity vector, $e _ { B }$ is the bottom erosion/ deposition rate, and $e _ { s , b - s }$ is the sediment mass exchange between bed and suspended load. The second-order tensor $\mathbf { U } _ { l } \mathbf { U } _ { l }$ (resp. $\mathbf { U } _ { s } \mathbf { U } _ { s } )$ represents the diadic product of the phaseaveraged water (resp. solid) velocity with itself. Finally, denoting with D the stress due to drag exchanged between the two phases, the source terms of momentum equations $\mathbf { S } _ { l }$ and $\mathbf { S } _ { s , b }$ are 

$$
\mathbf {S} _ {l} = \frac {\pmb {\tau} _ {B , l}}{\rho_ {l}} + \frac {\mathbf {D}}{\rho_ {l}}\tag{7}
$$

$$
\mathbf {S} _ {s, b} = \frac {\boldsymbol {\tau} _ {B , s}}{\rho_ {s}} - \frac {\mathbf {D}}{\rho_ {s}}\tag{8}
$$

in which $\tau _ { B , l }$ and $\pmb { \tau } _ { B , s }$ are the bottom shear stresses on the liquid and the solid phases, respectively. The drag force of the water on the solid particles, D, is evaluated as 

$$
\mathbf {D} = \rho_ {l} C _ {D} \frac {\delta_ {s , b}}{d} (\mathbf {U} _ {l} - \mathbf {U} _ {s}) | \mathbf {U} _ {l} - \mathbf {U} _ {s} |\tag{9}
$$

where $C _ { D } =$ bulk drag coefficient. The shear stress acting on the solid phase $\tau _ { B , s }$ is expressed as 

$$
\frac {\boldsymbol {\tau} _ {B , s}}{\rho_ {s}} = \mu_ {d} g \delta_ {s, b} \frac {r}{r + 1} \frac {\mathbf {U} _ {s}}{| \mathbf {U} _ {s} |} + \alpha \mathbf {U} _ {s} | \mathbf {U} _ {s} |\tag{10}
$$

in which $\mu _ { d }$ is the dynamic friction coefficient. Eq. (10) accounts for both frictional, expressed through Mohr-Coulomb law, and interparticle collisional (Bagnold 1956) stresses. Following Seminara et al. (2002), the shear stress on the liquid phase is evaluated by the following relation: 

$$
\boldsymbol {\tau} _ {B, l} = \rho_ {l} \frac {\mathbf {U} _ {l}}{C _ {\mathrm{Ch}} ^ {2}} | \mathbf {U} _ {l} | - \boldsymbol {\tau} _ {B, s} + \rho_ {l} g \delta_ {s, b} \mathbf {s} _ {B}\tag{11}
$$

where ${ \bf s } _ { B } = { \bf \frac { \pi } { \pi } }$ bottom slope. The first term is evaluated by means of the Chezy uniform flow formula, $C _ { C h }$ being the dimensionless Chezy coefficient. 

The bottom entrainment/deposition is expressed through the following formula proposed by Pontillo et al. (2010): 

$$
e _ {B} = w _ {s} \frac {T ^ {3 / 2} - C _ {s , b}}{1 - p}\tag{12}
$$

in which $w _ { s }$ denotes the sediment settling velocity and $C _ { s , b }$ is the bed load concentration. The dimensionless mobility parameter $T$ accounts for the excess of the mobilizing stresses onto the bottom surface with respect to the resisting ones (van Rijn 1984). A large number of experiments have shown that the settling velocity reduces as the particle concentration increases. The following semiempirical formula (Richardson and Zaki 1954) is therefore considered to evaluate the sediment settling velocity: 

$$
w _ {s} = w _ {t} (1 - C _ {s, b}) ^ {n}\tag{13}
$$

in which $w _ { t } =$ terminal settling velocity of a single particle in an indefinite fluid. According to Baldock et al. (2004), the exponent n is about 2.5 for particles with diameter of 1 mm, whereas it increases up to 5 for smaller sediments. 

The mobility parameter T is herein defined as 

$$
T = \frac {\left| \boldsymbol {\tau} _ {B , l} + \boldsymbol {\tau} _ {B , s} - \boldsymbol {\tau} _ {c} - \boldsymbol {\tau} _ {B} \right|}{\left| \boldsymbol {\tau} _ {c} + \boldsymbol {\tau} _ {B} \right|}\tag{14}
$$

where $\pmb { \tau } _ { c } =$ threshold shear stress for sediment motion and $| \pmb { \tau } _ { B } | =$ $\mu _ { s } r g \delta _ { s , b }$ is the Mohr-Coulomb stress at the bottom, with $\mu _ { s }$ the static friction coefficient. Under clear-water conditions, Eq. (12) states that the erosion rate scales with the $3 / 2$ power of the van Rijn transport parameter, which is consistent with van Rijn findings (van Rijn 1984). 

The solid exchange between the bed and suspended load is modeled through a first-order kinetic law (Wu et al. 2000) 

$$
e _ {s, b - s} = \beta \omega (C _ {s, s} ^ {*} - C _ {s, s})\tag{15}
$$

in which $C _ { s , s }$ represents the depth-averaged suspended sediment concentration, $C _ { s , s } ^ { * }$ is the corresponding capacity value. The exchange is modulated by $\beta$ and ω coefficients: the former relates the depth-averaged values to the local ones; the latter expresses the adaptation of suspended load and it is usually assumed as the sediment settling velocity $( \mathrm { i } . \mathrm { e } . , \omega = w _ { s } ) ,$ as it is also done herein. The expression proposed by Armanini and Di Silvio (1988) is used to evaluate $\beta .$ 

The capacity value for suspended sediment concentration is estimated through the following formula proposed by Wu et al. (2000) and Wu (2007): 

$$
C _ {s, s} ^ {*} = 0. 0 0 0 0 2 6 2 \frac {C _ {s , b} \sqrt {g d r} d}{| \mathbf {U} _ {l} | (C _ {s , b} h - \delta_ {s , b})} \left[ \left(\frac {\theta_ {0}}{\theta_ {c}} - 1\right) \frac {| \mathbf {U} _ {l} |}{w _ {s}} \right] ^ {1. 7 4}\tag{16}
$$

where $\theta _ { 0 } = \tau _ { 0 } / ( \rho _ { l } g d r )$ is the Shields parameter computed through the modulus of the shear stress $\tau _ { 0 }$ at the bed without considering the transport layer, and $\theta _ { c }$ is the corresponding threshold value for the sediment transport initiation. 

## Model Closures

The α and $C _ { D }$ coefficients may be estimated from existing empirical formulas (e.g., Maude and Whitmore 1958), which however introduce other parameters. As an alternative, in the present paper both coefficients are evaluated based on the analysis of uniform flow conditions. To this aim, the model is first applied to a uniform flow characterized by a bottom slope $s _ { B } .$ In such a condition, the two-phase conservation Eqs. (1)–(6) reduce to the following set of relations: 

$$
g (\delta_ {l} + \delta_ {s, b} + \delta_ {s, s}) s _ {B} = \frac {\tau_ {B , l}}{\rho_ {l}} + \frac {D}{\rho_ {l}}\tag{17}
$$

$$
g \delta_ {s, b} r s _ {B} = \frac {\tau_ {B , s}}{\rho_ {l}} - \frac {D}{\rho_ {l}}\tag{18}
$$

$$
C _ {s, b} = T ^ {3 / 2}\tag{19}
$$

$$
\beta C _ {s, s} = C _ {s, b}\tag{20}
$$

Similarly to Parker et al. (2003), the following scaling law for the bed load volume for unit bottom area is assumed: 

$$
\frac {\delta_ {s , b}}{d} = k _ {1} (\theta_ {0} - \theta_ {c})\tag{21}
$$

with $k _ { 1 }$ a dimensionless coefficient. Although Eq. (21) was deduced only for the low Shields parameter, i.e., $\theta _ { 0 } \leq 0 . 1$ (Fernandez-Luque and Van Beek 1976), recent experiments (Lajeunesse et al. 2010) have confirmed its validity up to $\theta _ { 0 }$ 0.2. In the present analysis, Eq. (21) is therefore applied even for a higher Shields number. 

The peculiarities of the solid particles motion in the bed load, through saltation, rolling, and sliding have been thoroughly investigated through experimental studies, which have suggested that sediment velocity is different from that of the carrying fluid. Several formulas have been proposed for its evaluation, witnessing the importance of its correct computation for bed load modeling. In particular, Meland and Norrman (1966) deduced an empirical expression of the sediment average transport velocity in terms of shear velocity, roughness size, and particle diameter based on a series of experiments with glass beads rolling on a bed of homogenously sized particles. The dimensional nature of this formula limits its validity to the range of the investigated experimental conditions. Fernandez-Luque and van Beek (1976), starting from experiments carried out with a loose bed, proposed the following expression of the particles average transport velocity $U _ { p } { \mathrm { : } }$ 

$$
U _ {p} = c _ {a} (u _ {*} - 0. 7 u _ {* c})\tag{22}
$$

in which $u _ { * } =$ shear velocity; $u _ { * c } =$ corresponding value in the Shields critical condition; and $\begin{array} { r l } { c _ { a } } & { { } = } \end{array}$ dimensionless constant approximately equal to 11.5. 

A theoretical consideration about the dynamics of the bed load sediment transport led Bridge and Dominic (1984) to deduce the following expression for the bed grain velocity: 

$$
U _ {g} = c _ {b} (u _ {*} - u _ {* c})\tag{23}
$$

with $c _ { b } = \sqrt { \tan \mu _ { d } } w _ { s } / u _ { * c } .$ 

Moreover, Sekine and Kikkawa (1992), presenting a deterministicprobabilistic model to investigate the nature of the bed load motion, proposed the following expression for the bed load layer averaged mean velocity of saltation: 

$$
\frac {U _ {m}}{\sqrt {g d r}} = 8 \frac {u _ {*}}{w _ {s}} \left(1 - \frac {u _ {* c}}{u _ {*}}\right) ^ {1 / 2}\tag{24}
$$

The effectiveness of the dimensionless parameters of Sekine and Kikkawa (1992) for describing the motion of sediment particles over transitionally rough beds has been successively confirmed by Papanicolaou et al. (2002b) and Ramesh et al. (2011). 

Seminara et al. (2002), in deriving an entrainment-based model of sediment transport that neither satisfies nor suffers from the drawbacks of the Bagnold constraint, proposed a slight modification of the Fernandez Luque and van Beek (1976) formula, which reads 

$$
U _ {p} = c _ {a} ^ {\prime} (\tau - \tau_ {c}) ^ {1 / 2}\tag{25}
$$

with the dimensionless coefficient $c _ { a } ^ { \prime }$ ranging between 8 and 9. Recently Julien and Bounvilay (2013), based on a dimensional and regression analysis carried out considering bed load particles on smooth and rough rigid plane surfaces, proposed a simple singleparameter relation, which expresses the bed load particle velocity in terms of the shear velocity and of the logarithm of the Shields parameter of the boundary roughness. 

In what follows, following Seminara et al. (2002), the solidphase average velocity in the bed load layer is assumed to be 

$$
\frac {U _ {s}}{\sqrt {g d r}} = k _ {2} (\theta_ {0} - \theta_ {c}) ^ {1 / 2}\tag{26}
$$

with $k _ { 2 }$ an experimental dimensionless coefficient. By postulating the validity of Eqs. (21) and (26), the following expression of the bed load solid discharge is deduced: 

$$
\frac {U _ {s} \delta_ {s , b}}{d \sqrt {g d r}} = k _ {1} k _ {2} \left(\frac {\tau_ {0} - \tau_ {c}}{\rho_ {l} g d r}\right) ^ {3 / 2} = k _ {1} k _ {2} (\theta_ {0} - \theta_ {c}) ^ {3 / 2}\tag{27}
$$

Eq. (27) has the same structure of the well-known Meyer-Peter and Müller (1948) formula, which is exactly reproduced provided that the $k _ { 1 } k _ { 2 }$ product is set equal to the Meyer-Peter and Müller coefficient $( K _ { \mathrm { M P M } } )$ $K _ { \mathrm { M P M } }$ ranges from about 4, as indicated in the reanalysis of original Meyer-Peter and Müller’s dataset described in Wong and Parker (2006), to 12, used in the numerical simulations reported in El Kadi Abderrezzak and Paquier (2011). The original and most used value of 8 (Meyer-Peter and Müller 1948) is adopted in what follows. Assuming the classical value $K _ { \mathrm { M P M } } = 8 ,$ , the two empirical parameters $k _ { 1 }$ and $k _ { 2 }$ are fixed by considering the bounds deriving by the consistency of the model, as shown in the following. 

The water velocity may be computed through Chezy’s law 

$$
\frac {U _ {l}}{\sqrt {g d r}} = C _ {C h} \theta_ {0} ^ {1 / 2}\tag{28}
$$

Finally, it is postulated that the shear stress acting on the liquid phase may be represented as follows: 

$$
\tau_ {B, l} = \tau_ {c} + c _ {1} (\tau_ {0} - \tau_ {c})\tag{29}
$$

with $c _ { 1 }$ a nonnegative parameter smaller than unity, i.e., $0 \leq c _ { 1 } \leq 1$ In fact, the case $c _ { 1 } = 0$ corresponds to the Bagnold’s hypothesis, i.e., the shear between fluid and bottom reduces to the critical value (Bagnold 1956). On the other hand, the condition $c _ { 1 } = 1$ implies that the shear stress acting on the liquid phase equals the corresponding value in absence of sediment transport, i.e., no momentum is transferred to the solid phase. However—as it will be shown later—a more restrictive upper bound may be specified for it. While clear indications may be found in the literature for estimating the $C _ { C h }$ and $K _ { \mathrm { M P M } }$ coefficients in their well-defined variability ranges, the dimensionless nonnegative coefficient $c _ { 1 }$ represents a free model parameter. In the “Results” section, classical literature values are assumed for $C _ { C h }$ and $K _ { \mathrm { M P M } } .$ , while the $c _ { 1 }$ coefficient is allowed to vary in order to investigate its influence on the model predictions. 

Substituting the relations (21), (26), (28), and (29) into $\operatorname { E q . }$ (17), the following expression of the drag coefficient may be easily obtained: 

$$
C _ {D} = \frac {1 - c _ {1}}{k _ {1}} \frac {\rho_ {l} g d r}{\left[ C _ {C h} \tau_ {0} ^ {1 / 2} - k _ {2} (\tau_ {0} - \tau_ {c}) ^ {1 / 2} \right] ^ {2}}\tag{30}
$$

The substitution of Eqs. (21) and (26) into the momentum equation of the solid phase in the bed load layer, Eq. (18), gives the following expression for α: 

$$
\alpha = \frac {(1 - c _ {1}) - k _ {1} (\mu_ {d} - s _ {B})}{(r + 1) k _ {2} ^ {2}}\tag{31}
$$

Expressions (30) and (31), strictly valid only in uniform flow, are used herein also in nonuniform conditions considering the local and instantaneous values of $s _ { B }$ and $\tau _ { 0 }$ for a fixed value of $c _ { 1 }$ . As far as the $c _ { 1 }$ value is concerned, Eq. (31) shows that the positivity of the α coefficient imposes the following upper bound: 

$$
c _ {1} \leq 1 - k _ {1} (\mu_ {d} - s _ {B})\tag{32}
$$

The considered closures suggest a way to select the value for the $k _ { 1 }$ coefficient, which has been experimentally found to vary between 0.66 (Seminara et al. 2002) and 2.51 (Lajeunesse et al. 2010). Indeed, rewriting the transport stage parameter T as 

$$
T = \frac {(\theta_ {0} - \theta_ {c}) (1 - k _ {1} \mu_ {s})}{\theta_ {c} + k _ {1} \mu_ {s} (\theta_ {0} - \theta_ {c})}\tag{33}
$$

and the concentration $C _ { s , b }$ as 

$$
C _ {s, b} = \frac {\delta_ {s , b}}{K _ {s} d}\tag{34}
$$

with $K _ { s }$ as the ratio of the bed load layer thickness to sediment diameter, the bottom entrainment/deposition condition (19) leads to the following expression for $K _ { s }$ : 

$$
K _ {s} = \frac {k _ {1}}{(1 - k _ {1} \mu_ {s}) ^ {3 / 2}} \frac {[ \theta_ {c} + k _ {1} \mu_ {s} (\theta_ {0} - \theta_ {c}) ] ^ {3 / 2}}{(\theta_ {0} - \theta_ {c}) ^ {1 / 2}}\tag{35}
$$

Moreover, accounting for Eqs. (21) and (35) may be equivalently rewritten in terms of the bed load volume for unit bottom area as follows: 

$$
K _ {s} = \frac {k _ {1} ^ {3 / 2}}{(1 - k _ {1} \mu_ {s}) ^ {3 / 2}} \frac {[ \theta_ {c} + \mu_ {s} \delta_ {s , b} / d ] ^ {3 / 2}}{(\delta_ {s , b} / d) ^ {1 / 2}}\tag{36}
$$

Eqs. (35) or (36) indicates that the positiveness of $K _ { s }$ implies the following condition on $k _ { 1 } \mathrm { { : } }$ : 

$$
k _ {1} <   \frac {1}{\mu_ {s}}\tag{37}
$$

Furthermore, for sufficiently large values of the shear stress, i.e., $( \theta _ { 0 } - \theta _ { c } ) \gg \theta _ { c }$ , as those corresponding to sheet-flow regime, Eq. (35) can be approximated as 

$$
K _ {s} ^ {S F} \cong \frac {k _ {1} ^ {5 / 2} \mu_ {s} ^ {3 / 2}}{(1 - k _ {1} \mu_ {s}) ^ {3 / 2}} (\theta_ {0} - \theta_ {c})\tag{38}
$$

and therefore the bed load concentration asymptotically approaches the value 

$$
C _ {s, b} ^ {S F} = \frac {(1 - k _ {1} \mu_ {s}) ^ {3 / 2}}{k _ {1} ^ {3 / 2} \mu_ {s} ^ {3 / 2}}\tag{39}
$$

Since the asymptotic concentration $\operatorname { E q . }$ (39) cannot exceed the sediment concentration in the erodible bottom, an additional condition for the $k _ { 1 }$ value has to be respected 

$$
k _ {1} \geq \frac {1}{\mu_ {s} [ 1 + (1 - p) ^ {2 / 3} ]}\tag{40}
$$

In what follows, the value of $k _ { 1 }$ is evaluated as the average between the lower Eq. (37) and upper Eq. (40) bounds 

$$
k _ {1} = \frac {1}{2 \mu_ {s}} \frac {2 + (1 - p) ^ {2 / 3}}{1 + (1 - p) ^ {2 / 3}}\tag{41}
$$

It is easy to verify that for common values of the porosity $( p )$ and of the static friction coefficient $( \mu _ { s } ) ,$ , Eq. (41) provides values for the $k _ { 1 }$ coefficient within the range of empirical values mentioned above. Furthermore, assuming the validity of the Meyer-Peter and Müller formula, the $k _ { 2 }$ coefficient is determined as 

$$
k _ {2} = \frac {K _ {\mathrm{MPM}}}{k _ {1}}\tag{42}
$$

In Fig. 1, the consistency of the above set of closures is verified by comparing the prediction of the dimensionless saltation height provided by Eq. (35), with available experimental (Lee and Hsu 1994; Nino et al. 1994; Nino and Garcia 1998; Lee et al. 2000) and numerical (Wiberg and Smith 1985) results. Since unfortunately the considered references do not specify the values of porosity and of the static friction coefficient, Eqs. (35) and (41) have been applied considering two reasonable pairs of $( \mu _ { s } , p ) .$ , namely, (0.5, 0.6) and (1.0, 0.4). On the other hand, accordingly with the values provided for the dimensionless threshold shear stress in the reference data, $\theta _ { c }$ has been assumed equal to 0.03 [Fig. 1(a)] in the comparison with data of Lee and Hsu (1994) and Wiberg and Smith (1985), and equal to 0.06 in the comparison with data of Lee et al. (2000), Nino and Garcia (1998) and Nino et al. (1994), [Fig. 1(b)]. 

Fig. 1 shows that Eqs. (35) and (41) provide relatively accurate predictions of the bed load layer thickness up to values of the Shields parameter order of unity. The fairly good agreement justifies the use of the relation (21) for the sediment volume for unit bottom area in combination with the entrainment formulation proposed by Pontillo et al. (2010) up to $\theta _ { 0 } \approx 1$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/76f2b2dcdca31711eaa0f353ca33643ca855220da8d80124da386e394b6f545d.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/b5021371057e70ce121b7a2e7e3d7edfabe059dbe14d8bc93059b3addcb899b7.jpg)



Fig. 1. Comparison between predictions by Eq. (35) and literature data: (a) $\theta _ { c } = 0 . 0 3 ;$ (b) $\theta _ { c } = 0 . 0 6$


## Model Properties and Numerical Method

In order to show the hyperbolic character of the presented flow model, system Eqs. (1)–(6) is rewritten in quasi-linear form. Accounting for Eqs. (34) and (36) and without considering the source terms, it reads 

$$
\mathbf {C} \frac {\partial \mathbf {W}}{\partial t} + \mathbf {A} \frac {\partial \mathbf {W}}{\partial x} + \mathbf {B} \frac {\partial \mathbf {W}}{\partial y} = 0
$$

<sub>ð</sub><sup>43</sup><sub>Þ</sub> 

in which, denoting with U and V the x and y components of velocity vector for both phases, the unknowns’ vector W is 

$$
\mathbf {W} = \left[ \begin{array}{l} \delta_ {l} \\ U _ {l} \\ V _ {l} \\ \delta_ {s, b} \\ U _ {s} \\ V _ {s} \\ z _ {B} \\ \delta_ {s, s} \end{array} \right]\tag{44}
$$

and the C, A, B, matrices may be easily deduced from Eqs. (1)–(6), through standard algebra. 

Following Courant and Hilbert (1961), the mathematical character of system (43) is investigated by looking for the eigenvalues of the matrix 

$$
\mathbf {M} = \mathbf {C} ^ {- 1} (\mathbf {A} n _ {x} + \mathbf {B} n _ {y})\tag{45}
$$

with $n _ { x }$ and $n _ { y }$ the director cosines of an arbitrary direction in the (x; y) plane of the unitary vector n. The eigenvalues read 

$$
\begin{array}{l} \lambda_ {1} = 0 \quad \lambda_ {2, 3} = \mathbf {U} _ {l} \cdot \mathbf {n} \quad \lambda_ {4} = \mathbf {U} _ {s} \cdot \mathbf {n} \\ \lambda_ {5, 6} = \mathbf {U} _ {l} \cdot \mathbf {n} \pm \sqrt {g h} \sqrt {\frac {\delta_ {l} + \delta_ {s , s}}{\delta_ {l}}} \\ \lambda_ {7, 8} = \mathbf {U} _ {s} \cdot \mathbf {n} \pm \sqrt {\frac {g d r}{2 (r + 1)}} \sqrt {K _ {s} + \frac {d K _ {s}}{d \delta_ {s , b}} \delta_ {s , b}} \end{array}\tag{46}
$$

in which the derivative of the dimensionless bed load layer thickness with respect to $\delta _ { s , b }$ has the following expression: 

$$
\frac {d K _ {s}}{d \delta_ {s , b}} = \frac {1}{2 \delta_ {s , b}} \left(\frac {k _ {1}}{1 - k _ {1}}\right) ^ {3 / 2} \sqrt {\frac {d}{\delta_ {s , b}} \theta_ {c} + \mu_ {s}} \left(2 \mu_ {s} \frac {\delta_ {s , b}}{d} - \theta_ {c}\right)\tag{47}
$$

Accounting for (47) eigenvalues $\lambda _ { 7 , 8 }$ may be equivalently rewritten as follows: 

$$
\begin{array}{l} \lambda_ {7, 8} = \mathbf {U} _ {s} \cdot \mathbf {n} \pm \frac {1}{2} \sqrt {\frac {g d r}{r + 1}} \\ \times \sqrt {\left(\frac {k _ {1}}{1 - k _ {1}}\right) ^ {3 / 2} \left(4 \mu_ {s} \frac {\delta_ {s , b}}{d} + \theta_ {c}\right) \sqrt {\frac {d}{\delta_ {s , b}} \theta_ {c} + \mu_ {s}}} \end{array}\tag{48}
$$

From Eqs. (46) and (48) it follows that, independently of the n unitary vector, the matrix M possesses only real eigenvalues. Therefore, the present two-phase model is always hyperbolic, and the characteristics theory allows to define the correct number of conditions on each boundary of the computational domain. 

The model represented by Eqs. (1)–(6) may be equivalently rewritten in a compact form as follows: 

$$
\frac {\partial \mathbf {U} _ {c}}{\partial t} + \frac {\partial \mathbf {F} (\mathbf {U} _ {c})}{\partial x} + \frac {\partial \mathbf {G} (\mathbf {U} _ {c})}{\partial y} + \mathbf {N} + \mathbf {S} _ {c} = 0\tag{49}
$$

in which 

$$
\begin{array}{l} \mathbf {U} _ {c} = \left( \begin{array}{c} \delta_ {l} \\ \delta_ {s, b} \\ \delta_ {s, s} \\ U _ {l} \delta_ {l} \\ V _ {l} \delta_ {l} \\ U _ {s} \delta_ {s, b} \\ V _ {s} \delta_ {s, b} \\ z _ {B} \end{array} \right); \qquad \mathbf {N} = \left( \begin{array}{c} 0 \\ 0 \\ 0 \\ g (\delta_ {l} + \delta_ {s, b} + \delta_ {s, s}) \frac {\partial z _ {B}}{\partial x} \\ g (\delta_ {l} + \delta_ {s, b} + \delta_ {s, s}) \frac {\partial z _ {B}}{\partial y} \\ g \frac {r}{r + 1} \delta_ {s, b} \frac {\partial z _ {B}}{\partial x} \\ g \frac {r}{r + 1} \delta_ {s, b} \frac {\partial z _ {B}}{\partial y} \\ 0 \end{array} \right); \\ \mathbf {S} _ {c} = \left( \begin{array}{c} p e _ {B} \\ (1 - p) e _ {B} - e _ {s, b - s} \\ e _ {s, b - s} \\ S _ {l, x} \\ S _ {l, y} \\ S _ {s, x} \\ S _ {s, y} \\ e _ {B} \end{array} \right) \end{array}\tag{50}
$$

and 

$$
\begin{array}{l} \mathbf {F} = \left( \begin{array}{c} \delta_ {l} U _ {l} \\ \delta_ {s, b} U _ {s} \\ \delta_ {s, s} U _ {l} \\ \delta_ {l} U _ {l} ^ {2} + g \frac {\left(\delta_ {l} + \delta_ {s , b} + \delta_ {s , s}\right) ^ {2}}{2} \\ \delta_ {l} U _ {l} V _ {l} \\ \delta_ {s, b} U _ {s} ^ {2} + g \frac {r}{r + 1} \frac {\delta_ {s , b} ^ {2}}{2 C _ {s , b}} \\ \delta_ {s, b} U _ {s} V _ {s} \\ 0 \end{array} \right) \\ \mathbf {G} = \left( \begin{array}{c} \delta_ {l} V _ {l} \\ \delta_ {s, b} V _ {s} \\ \delta_ {s, s} V _ {l} \\ \delta_ {l} U _ {l} V _ {l} \\ \delta_ {l} V _ {l} ^ {2} + g \frac {\left(\delta_ {l} + \delta_ {s , b} + \delta_ {s , s}\right) ^ {2}}{2} \\ \delta_ {s, b} U _ {s} V _ {s} \\ \delta_ {s, b} V _ {s} ^ {2} + g \frac {r}{r + 1} \frac {\delta_ {s , b} ^ {2}}{2 C _ {s , b}} \\ 0 \end{array} \right) \end{array}\tag{51}
$$

Vector N represents the nonconservative terms in the partial differential system, arising from the bed slope source term. 

The system (49) can be solved with any of the numerical schemes commonly used for hyperbolic partial differential equations (PDEs). The finite volume solver used in (Leopardi et al. 2002; Greco et al. 2012a) has been adapted to solve the PDEs of the two-phase model, along with an appropriate treatment of the bed slope source term N (Valiani and Begnudelli 2006; Greco et al. 2008a). To this aim, with reference to a structured rectangular mesh Eq. (49) is written in the following semidiscrete conservative form 

$$
\frac {d \bar {\mathbf {U}} _ {c}}{d t} = - \frac {1}{A _ {0}} \left[ \sum_ {k = 1} ^ {4} \left(\mathbf {H} _ {k} \cdot l _ {k} \mathbf {n} _ {k}\right) - \bar {\mathbf {S}} _ {c} \right]\tag{52}
$$

In Eq. (52), the overbar denotes the averaging over the computational cell of area $A _ { 0 } , l _ { k }$ is the length of the kth side of the cell, ${ \bf n } _ { k }$ is the normal vector and $\mathbf { H } _ { k }$ is the average value of the flux on the same side, defined as 

$$
\mathbf {H} _ {k} = \mathbf {F} ^ {\prime} n _ {x} + \mathbf {G} ^ {\prime} n _ {y}\tag{53}
$$

being $\mathbf { F } ^ { \prime }$ and $\mathbf { G } ^ { \prime }$ the vectors of the numerical fluxes, modified as follows to include the slope terms: 

$$
\begin{array}{l} \mathbf {F} ^ {\prime} = \left( \begin{array}{c} \delta_ {l} U _ {l} \\ \delta_ {s, b} U _ {s} \\ \delta_ {s, s} U _ {l} \\ \delta_ {l} U _ {l} ^ {2} + g \frac {(\delta_ {l} + \delta_ {s , b} + \delta_ {s , s})}{2} [ (\delta_ {l} + \delta_ {s, b} + \delta_ {s, s}) + z _ {B} - \tilde {z} ] \\ \delta_ {l} U _ {l} V _ {l} \\ \delta_ {s, b} U _ {s} ^ {2} + g \frac {r}{r + 1} \frac {\delta_ {s , b}}{2 C _ {s , b}} [ \delta_ {s, b} + z _ {B} - \tilde {z} ] \\ \delta_ {s, b} U _ {s} V _ {s} \\ 0 \end{array} \right) \\ \mathbf {G} ^ {\prime} = \left( \begin{array}{c} \delta_ {l} V _ {l} \\ \delta_ {s, b} V _ {s} \\ \delta_ {s, s} V _ {l} \\ \delta_ {l} U _ {l} V _ {l} \\ \delta_ {l} V _ {l} ^ {2} + g \frac {(\delta_ {l} + \delta_ {s , b} + \delta_ {s , s})}{2} [ (\delta_ {l} + \delta_ {s, b} + \delta_ {s, s}) + z _ {B} - \tilde {z} ] \\ \delta_ {s, b} U _ {s} V _ {s} \\ \delta_ {s, b} V _ {s} ^ {2} + g \frac {r}{r + 1} \frac {\delta_ {s , b}}{2 C _ {s , b}} [ \delta_ {s, b} + z _ {B} - \tilde {z} ] \\ 0 \end{array} \right) \end{array}\tag{54}
$$

z~ is the bed elevation at the side of the cell opposite the one on which flux has to be evaluated; the terms in the square bracket are considered null if negative (Greco et al. 2008a). 

Time integration of Eq. (52) is performed with a predictorcorrector (McCormack) scheme 

$$
\begin{array}{l} \bar {\mathbf {U}} _ {c} ^ {*} = \bar {\mathbf {U}} _ {c} ^ {t} - \frac {\Delta t}{A _ {0}} \left[ \sum_ {k = 1} ^ {4} (\mathbf {H} _ {k} ^ {t} \cdot l _ {k} \mathbf {n} _ {k}) - \bar {\mathbf {S}} ^ {t} \right] \\ \bar {\mathbf {U}} _ {c} ^ {* *} = \bar {\mathbf {U}} _ {c} ^ {t} - \frac {\Delta t}{A _ {0}} \left[ \sum_ {k = 1} ^ {4} (\mathbf {H} _ {k} ^ {*} \cdot l _ {k} \mathbf {n} _ {k}) - \bar {\mathbf {S}} ^ {*} \right] \end{array} \quad \bar {\mathbf {U}} _ {c} ^ {t + \Delta t} = \frac {\bar {\mathbf {U}} _ {c} ^ {*} + \bar {\mathbf {U}} _ {c} ^ {* *}}{2}\tag{55}
$$

The numerical fluxes at the interfaces are computed by a threepoint parabolic interpolation of the conserved variables values. In the predictor stage, two cells on a side of the interface and one on the opposite side are considered, vice versa in the corrector stage. The numerical stability of the proposed method is guaranteed provided that the Courant–Friedrichs–Lewy condition is satisfied for the largest eigenvalue [Eq. (46)]. 

## Test Cases and Results

In the next two subsections the proposed model is tested against two laboratory experiments: a one-dimensional dam break, over a dry erodible bed (Capart and Young 1998), and a two-dimensional dam break, over both dry and wet bed (Soares-Frazão et al. 2012). Finally, in the last section of this paragraph, the present model is compared to four existing non-equilibrium models. 

## One-Dimensional Dam Break

The first test case is the fast geomorphic transient experimentally investigated by Capart and Young (1998). The experiments were carried out at National Taiwan University, and they consist of small-scale laboratory dam break of initial water depth $h _ { 0 } = 1 0$ cm over an erodible bed in a prismatic rectangular channel. Notably, a very light sediment was used (density $\rho _ { s } = 1 , 0 4 8 \mathrm { ~ k g } \mathrm { m } ^ { - 3 } )$ with $d = 6 . 1$ mm. Scouring propagates both upstream and downstream of the dam, where intense erosion occurs. Apart from the near-field evolution soon after the dam removal, the flood wave exhibits a rather regular shape characterized by a steep sediment-laden bore, at the front of the wave, and an enduring weak hydraulic jump at the center of the wave. 

As indicated by the experimenters, the bottom porosity p has been fixed equal to $0 . 6 ,$ , whereas the sediment free-fall velocity $w _ { t }$ in Eq. (13) is assumed equal to 0.067 $\mathrm { m } / \mathrm { s }$ . The settling velocity $w _ { s }$ is computed through Eq. (13) at each point and time accordingly to the actual concentration value and with the n value fixed equal to 2.5. The values of the static and dynamic friction coefficients are $\mu _ { s } = 0 . 5 2$ and $\mu _ { d } = 0 . 3 2$ , respectively. The dimensionless Chezy coefficient has been evaluated by Griffiths’ (1981) formula for a value of the $h / d$ ratio of about 12. The threshold Shields number was fixed at the classical value of $\theta _ { c } = 0 . 0 4 7$ and the Meyer-Peter and Müller coefficient $( K _ { \mathrm { M P M } } )$ has been assumed equal to 8. The $k _ { 1 }$ and $k _ { 2 }$ coefficients have been evaluated through Eqs. (41) and (42), respectively, and their values are $k _ { 1 } = 1 . 0 5$ and $k _ { 2 } = 7 . 6 2$ . Finally, the upper bound value of the free parameter $c _ { 1 }$ , deduced by Eq. (32) is 0.44. 

Simulations have been carried out with a grid size $\Delta x =$ 0.010 m and $\Delta t = 1 / 4 , 0 9 6 \ \mathrm { s } .$ The computational domain was sufficiently long to exclude any influence of the boundary conditions. Three different values of the $c _ { 1 }$ parameter, namely $c _ { 1 } = 0 , c _ { 1 } = 0 . 2$ and $c _ { 1 } = 0 . 4$ , have been considered. In Fig. 2 two snapshots of the experimental results from Fraccarollo and Capart (2002), corresponding to $t = 0 . 4$ s and $t = 0 . 5 ~ \mathrm { s }$ after dam removal, are compared with the computed results. The numerical results show a very limited sensitivity to the $c _ { 1 }$ value and moreover they indicate that the model predictions closely agree with the main features of the process, i.e., the celerity of the downstream tail, the free surface profile upstream and downstream the dam, and the scour of the bottom. The shape of the scour strongly resembles the experimental one, with a steep adverse slope just downstream the original dam location $( x = 0 )$ , followed by a nearly horizontal scoured bed. A general slight underestimation of the maximum scour occurring just upstream the bore is however observed at $t = 0 . 4 ~ \mathrm { s }$ . The observed weak hydraulic jump is also qualitatively reproduced in the simulations, with a bore appearing more upstream than in the experiments and with a sharper front. 

As far as the sediment transport reproduction is concerned, Fig. 3(a) depicts in the space-time plane the suspended sediment discharge values $q _ { s , s } = \delta _ { s , s } U _ { l }$ divided by the total solid discharge $\begin{array} { r } { q _ { s , \mathrm { t o t } } = \delta _ { s , s } U _ { l } + \delta _ { s , b } U _ { s } , } \end{array}$ , while Fig. 3(b) reports the space-time evolution of the ratio $K _ { s } d / h$ . Even if in a large portion of the plane the suspended transport represents a small percentage, about 2%, of the total solid discharge, the map shows that there are some areas in which it increases up to 20%. The suspended solid discharge represents an appreciable contribution to the solid discharge only in a limited portion of the (x; t) plane, while it is absent in most of the region downstream to the original dam $( \mathrm { i } . \mathrm { e } . , x > 0 )$ , although in this region the Rouse number is less than 1 (results not shown herein). Such a result may be explained accounting for that, downstream the original dam position, the bed load thickness saturates the full flow depth [Fig. 3(b)] and therefore the solid discharge is entirely conveyed as bed load. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/8867be647334f413535a0b73ed8449c9f41c34d0458662e2a9451164c616c902.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/3c7c98555457ec50fc549b9da9125845cd29da30ecd357c663e9eda1c3483508.jpg)



Fig. 2. One-dimensional dry-bed test; measured and computed free-surface and bottom profiles: $\left( \mathrm { a } \right) t = 0 . 4 \ \mathrm { s } ; \left( \mathrm { b } \right) t = 0 . 5$ s after dam removal


## Two-Dimensional Dam Break

An example of a two-dimensional fast geomorphic transient involving a wide range of the Shields parameter values is provided by the experiments carried out within the NSF-Pire project (Soares-Frazão et al. 2012). 

The tests concern dam break waves expanding over a flat mobile bed, in a 3.6 m wide, 36 m long flume, whose geometry is reported in Fig. 4. The breached dam is represented by two impervious blocks and a 1.0 m wide gate located between the blocks. The sudden rise of the gate induces a flood wave expanding along both longitudinal and transversal directions. An initial 85-mm thick layer of coarse sand was put down upon the fixed bed, from 1 m upstream to 9 m downstream the gate. Sediments were constituted of a uniformly graded sand with $d = 1 . 6 1 1 0 ^ { - 3 }$ m with relative density $r = 1 . 6 3 ,$ with a bottom porosity $p = 0 . 4 2 .$ . The sediment free-fall velocity $w _ { t }$ is 0.18 $\mathrm { m } / \mathrm { s }$ . Also in this test case the settling velocity $w _ { s }$ has been computed through Eq. (13) with $n = 2 . 5$ and considering the actual concentration value. The following values of friction coefficients have been assumed $\mu _ { s } = 0 . 7 3$ and $\mu _ { d } = 0 . 6 3$ . The value of the $k _ { 1 }$ coefficient through Eq. (41) is $k _ { 1 } = 1 . 0 9$ . The threshold Shields parameter and the Meyer-Peter and Müller coefficients have been fixed equal to $\theta _ { c } = 0 . 0 4 7$ and $K _ { \mathrm { M P M } } = 8 .$ , as in the previous test, so that $k _ { 2 } = 7 . 3 4$ . The dimensionless Chezy coefficient has been similarly evaluated using the Griffiths’ formula. Here the ratio $h / d$ is about 200. The upper bound of the $c _ { 1 }$ parameter is 0.29. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/cbc66eae7310f1ec69af7ed1edacf10c0697a0a6bb388fee86ac91c46f399890.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/a38efbbb4439ed602652273fb0127fb1ed05e556d1233d5135a009c8491aaf40.jpg)



Fig. 3. Space-time maps of (a) suspended to total solid load ratio; (b) bed load thickness to flow depth ratio


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/175232aeec2c61ed4f2f2cedb8635f5240e4ebc359c0e043aed5943bc5fa0f58.jpg)



Fig. 4. NSF-PIRE Benchmark; scheme of the experimental setup (reprinted from Soares-Frazão et al. 2012, © International Association for Hydro-Environment Engineering and Research, reprinted by permission of Taylor and Francis Limited, http://www.tandfonline.com, on behalf of International Association for Hydro-Environment Engineering and Research)



Table 1. Gauges Locations for Test 1 Dry-Bed Test


<table><tr><td>Gauge n°</td><td>x (m)</td><td>y (m)</td></tr><tr><td>1</td><td>0.64</td><td>-0.5</td></tr><tr><td>2</td><td>0.64</td><td>-0.165</td></tr><tr><td>3</td><td>0.64</td><td>0.165</td></tr><tr><td>4</td><td>0.64</td><td>0.5</td></tr><tr><td>5</td><td>1.94</td><td>-0.99</td></tr><tr><td>6</td><td>1.94</td><td>-0.33</td></tr><tr><td>7</td><td>1.94</td><td>0.33</td></tr><tr><td>8</td><td>1.94</td><td>0.99</td></tr></table>


Table 2. Gauges Locations for Test 1 Wet-Bed Test


<table><tr><td>Gauge n°</td><td>x (m)</td><td>y (m)</td></tr><tr><td>1</td><td>0.64</td><td>-0.5</td></tr><tr><td>2</td><td>0.64</td><td>-0.165</td></tr><tr><td>3</td><td>0.64</td><td>0.165</td></tr><tr><td>4</td><td>0.64</td><td>0.5</td></tr><tr><td>5</td><td>2.34</td><td>-0.99</td></tr><tr><td>6</td><td>2.34</td><td>-0.33</td></tr><tr><td>7</td><td>2.34</td><td>0.33</td></tr><tr><td>8</td><td>2.34</td><td>0.99</td></tr></table>

Two configurations were experimentally investigated: (1) an initial water level of 47 cm in the upstream reservoir and no water downstream (dry-bed test); (2) an initial water level of 51 cm in the upstream reservoir and a water level of 15 cm downstream (wet-bed test). The time evolution of the water level was measured at eight gauges by means of ultrasonic probes (Fig. 4), whose location is indicated in Tables 1 and 2 for dry and wet bed test, respectively. The final topography was measured by a bottom profiler with 5 cm resolution along y. Further details about the experimental procedure may be found in the paper by Soares-Frazão et al. (2012). 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/ec8113df6b0b8db759436a72d555280bd3791002c9c92049a480865c48343957.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/e21d90d95cfc0dd16d5f739234167d743fa2f4b9ddb3d9dd14199395c951b873.jpg)


Both the dry-and wet-bed experiments have been simulated by means of a non-uniform mesh of about 35,000 cells, with variable size in x and $y$ directions. The smallest cells, used to discretize the erodible floodplain, have size $\Delta x = \Delta y = 2 . 5 \cdot 1 0 ^ { - 2 } \mathrm { ~ m ~ }$ . The adopted time step was $\Delta t = 1 / 2$ , 048 s. Freefall has been considered at the outlet section of the flume, whereas impervious boundaries have been considered for the flume sidewalls. 

With reference to Test Case 1 (dry-bed), Fig. 5 compares measured and computed time series of free-surface elevation at the gauge points, obtained with three different values of the $c _ { 1 }$ parameter, namely, 0, 0.1 and 0.2. Measures from symmetrical gauge points are grouped on the same plot. 

An estimate of the experiment reproducibility has been provided by Soares-Frazão et al. (2012) resulting in mean observed standard deviation between $\sigma _ { \mathrm { m e a n } } = 0 . 0 0 6 \div 0 . 0 1 6$ m with maximum values being between $\sigma _ { \operatorname* { m a x } } = 0 . 0 1 8 \div 0 . 0 3 2$ m, depending on the considered gauge. It is noticed that in all the gauges the arrival time of the surge caused by the dam failure is well captured, along with the general trend of the free-surface elevation decay after the surge transition. 

The experimental and simulated final bottom topographies for three values of the y coordinate $( y = 0 . 2 \ \mathrm { m } , \ y = 0 . 7$ m and $y = 1 . 4 5 \ \mathrm { m } )$ are compared in Fig. 6, still considering the same three different $c _ { 1 }$ values of Fig. 5. A slight but systematic underprediction of the deposition is observed in the simulated profile. This performance appears satisfactory if the scattering between the results of different repeated experimental runs is accounted for. Indeed, Soares-Frazão et al. (2012) estimated mean and maximum standard deviation of $\sigma _ { \mathrm { m e a n } } = 0 . 0 0 8$ m and $\sigma _ { \mathrm { m a x } } = 0 . 0 2 9$ m, respectively, with the latter value referring to the most intensely scoured zone. Moreover, the results depicted in both Figs. 5 and 6 confirm the limited influence of the $c _ { 1 }$ parameter on the results quality. 

Fig. 7 reports the vector plot of both water and sediment velocities at different times $( t = 2 \ \mathrm { s } , t = 5 \ \mathrm { s } , t = 2 0 \ \mathrm { s } ) .$ , showing the differences between the velocity fields of the two phases. In particular, the different alignment of the velocities vectors of the two phases is 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/8206aa376d98ecd732a1b567215aa7f5f89d75cd2d988f641adce67edcaf8d90.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/b9e03c76cb35cb4268fff64c5dbb15ba28f96643005544eebebfd17535ca3ead.jpg)



Fig. 5. Two-dimensional dry-bed test; measured and computed time series of free-surface elevation: (a) gauges 1 and $4 ;$ (b) gauges 2 and 3; (c) gauges 5 and 8; (d) gauges 6 and 7


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/f970e8b447496e93fe47941bb0b312ebf3dbabb5ff8b19024962a6b9c26116ee.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/e6214ece28b14fdf47fd3ec0f56b278691be00a86d179dc5d7da78d45ed74661.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/2039a0fb3f09f90ed0a957b4711b7c0116fbf8b669c23ca6e28589ff23891344.jpg)



Fig. 6. Two-dimensional dry-bed test; measured and simulated final bottom profiles: (a) $y = 0 . 2$ m; (b) $y = 0 . 7$ m; (c) y 1.45 m


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/d5ae39ffe37064eb84d981039b7364499880bbe85e0a4a684ef1554f2f1be96f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/0cc298e1cfbe718006bf1448b3e07892a0ff723cfe53a2e362e29f3792f88f54.jpg)



(c)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/641b0574ebfa2f1c62d38f065ba74cf4438206467e7b72d04522380fdeb7c00f.jpg)



x (m)



Fig. 7. Two-dimensional dry-bed test; velocity vector plot: (a) t  2 s; (b) t  5 s; (c) t  20 s after dam removal


evident for t 5 s, after that the flood wave impacted the sidewall and it was reflected toward the channel axis. The fluid flow is more responsive than the sediment to the impact of the wave. As far as the far-field t  20 s snapshot is considered, the sediment transport has ceased in the recirculation zone past the rigid blocks. Moreover, the symmetry of the velocity vectors respect to the longitudinal axis confirms the ability of the adopted numerical scheme to predict symmetric results. With reference to the same instants, the wide range of the Shields parameter of this flow is witnessed in Fig. 8. 

Finally, Fig. 9 represents the instantaneous values of $C _ { s , b }$ for the same times of Fig. 7. At all times, a steep transversal gradient of the concentration is observed in the narrow channel between the blocks. For $t = 2 \ \mathrm { s } ,$ the bulb-like flood wave exhibits a nearly constant concentration in its body and a gradual decrease close to the wave tip region, where the solid phase is transferred toward the suspension. However, maximum observed $C _ { s , s }$ values are smaller by more than one order of magnitude than the $C _ { s , b }$ ones (not reported). The results of both Figs. 8 and 9 also show a symmetric behavior respect to the longitudinal axis. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/2bea88257caabb10f34f3668fc95bc81538a0e784efe562f4008537515b4bbb4.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/8d756860d51ba381fab84e098e61919ef3f146940e4f07f533b1c3d7e1b8e0ac.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/439e859527da208b08a6f3465eedac9f3b70d351cfa524eeaf5a05a25266185f.jpg)



Fig. 8. Two-dimensional dry-bed test; Shields parameter distribution: (a) t  2 s; (b) t  5 s; (c) t  20 s after dam removal


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/56efd59bcf14c99ad2882e06a1a6ce5885b910164d2fa37a37a046bf1b800887.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/306a5d8b8641f7bd75a992f8255318565cab7ef1752ea314c1aab3f5e2adf293.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/1980a705aa184d4943165bd026055d52fc7c5e973d3e2592ba5088f41c8e22f3.jpg)



Fig. 9. Two-dimensional dry-bed test; bed load concentration distribution: (a) $t = 2 \ \mathrm { s } ;$ (b) t 5 s; (c) t 20 s after dam removal



t (s)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/7f5904eec011c12340f919d00e11dae92cd0a81e3cfd79741eb888a5f46cc78d.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/726cb29e0204bc72ca6c2e5ff2fe74bec135fed3ce63d37900a00e61f51e3ead.jpg)



(c)



(b)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/65ead0e991be3c2821ce5b9f41f6477f637b44ccb7484f8cce8c1fa8ae90378a.jpg)



t (s)



(d)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/ee7864775dbfd0f1dc1fc23b76f959b4f8b79a223f0a80d09b24a0efd8f62744.jpg)



Fig. 10. Two-dimensional wet-bed test; measured and computed time series of free-surface elevation: (a) gauges 1 and 4; (b) gauges 2 and 3; (c) gauges 5 and 8; (d) gauges 6 and 7


With reference to Test Case 2 (wet bed), Figs. 10 and 11 report the time series of the free-surface elevation at the different gauge points and of the final topography for the three longitudinal sections $y = 0 . 2 , 0 . 7$ and 1.45 m, respectively. The sensitivity respect to the $c _ { 1 }$ parameter is also represented. The results show that the present model is able to reproduce satisfactorily even in this test the wave propagation process (Fig. 10), independently of the $c _ { 1 }$ value. Moreover, the computed bed profile (Fig. 11) is characterized by bedforms in the scour hole with a comparable length than in the experiments, whereas the remaining of the profile is less wavy compared than the experimental one. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/99fea1d5be6891b7f9480881cbe8db9d3612db30c8ca3f30ff2061be4efa7f90.jpg)


The vector plot of both water and sediment velocities at different instants $( t = 2 \ \mathrm { s } , t = 5 \ \mathrm { s } , t = 2 0$ s) are represented in Fig. 12. As far as the direction of the liquid and solid velocity is concerned, the presence of the water downstream the dam tends to dampen the differences. On the other hand, the initial quiescent water downstream the dam obstacles the momentum diffusion, which leads to a significantly different shear stress distribution with respect to the dry-bed test case. Indeed, whereas the range of the shear stress values encountered by the flow is comparable with that of the previous test case, the spatial distribution is characterized by a more pronounced shear stress concentration in the region downstream the corner, as shown in Fig. 13. 

Along with the different shear stress distribution, the wet-bed test case differs significantly from the dry bed one also for the bed load concentration distribution. To enlighten such an aspect, the $C _ { s , b }$ distribution is represented in Fig. 14 with reference to the same instants considered for the previous case. At the first snapshot $\left( t = 2 ~ \mathrm { s } \right)$ , in fact, spatial gradients are more pronounced than in the dry-bed test case. At $t = 5 \ \mathrm { s } ,$ the $C _ { s , b }$ distribution is characterized by concentrations progressively reducing in the positive x direction. The nonuniform distribution evolves in time toward a more homogeneous one. In the near field, however, the capability of the present model to account for variable concentration seems fundamental for the bed load sediment routing. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/73776d2aa62d786be2afabe642d0162b42a02c1e7e28115a25fb734650e83461.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/4f0007bada3ee992174995ca94eae44b99ccc6633320f31f768f6b13f1488a5d.jpg)



Fig. 11. Two-dimensional wet-bed test; measured and simulated final bottom profiles: (a) $y = 0 . 2$ m; (b) $y = 0 . 7$ m; (c) y <sub>¼</sub> 1.45 m


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/f3f9a20fd488ab14fa8bb7116844743505d0b2b0e954fc73024cf71af4f3e807.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/23f8ab8e10a0b3a2390fb2e3203b2c5f2ab50570d8af89ad9a6fffc80ff76c61.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/d07454fa9c318cfe2f69feb0c6760da7f4c2fed87b75e119a87e476405aedc44.jpg)



Fig. 12. Two-dimensional wet-bed test; velocity vector plot: (a) $t = 2 \ \mathrm { s } ;$ (b) t <sub>¼</sub> 5 s; (c) t <sub>¼</sub> 20 s after dam removal


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/6ef5b3f7b0445767d74932c70931e3a091c46efa64c2f30022336e289bf64ca0.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/50f95897ee6c34ab4bb9bf9b2fdef1acba6f6e935a9685be6a101de4a29f8fc1.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/59603bf6290b5badb29e608f7ef45c93f8f6eea5e1753a7089277b69825b2bc2.jpg)



Fig. 13. Two-dimensional wet-bed test; Shields parameter distribution: (a) $t = 2 \ \mathrm { s } ;$ (b) t 5 s; (c) t 20 s after dam removal


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/3e310e3c90ac6ea3b11dabf7f2e20ff236899356199ad3bcd8fffa3b37c78819.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/b16c2f3e72259e48393bd9bf1a46176401f7da205855e09f5edd77de2f312173.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/acd9ea48732d40e09347817ecf4a775e80201969a35b6c8dad593702b31cbea8.jpg)



Fig. 14. Two-dimensional wet-bed test; bed load concentration distribution: (a) t <sub>¼</sub> 2 s; (b) t <sub>¼</sub> 5 s; (c) t <sub>¼</sub> 20 s after dam removal


## Comparison with Models in the Literature

In this section, the results of present model are compared with those obtained with four different models discussed in the literature review. 

The comparison concerns the main underlying assumptions of the different models, the evaluation of their specific parameters, the computational complexity (herein intended as the number of equations to be solved), along with the agreement with the experimental tests considered in the previous sections. 

As detailed in the “Model Closures” section, the present model essentially contains three dimensionless parameters, i.e. $C _ { C h } ,$ $K _ { \mathrm { M P M } }$ , and $c _ { 1 }$ . The parameters $C _ { C h }$ and $K _ { \mathrm { M P M } }$ may be evaluated based on extensive literature indications, whereas for $c _ { 1 }$ lower and upper bounds can be estimated. As far as the computational complexity is concerned, the one-dimensional (resp. two-dimensional) form of the proposed model needs the solution of five (resp. seven) differential equations expressing conservation principles of mass and momentum. Additionally, the bed evolution equation [Eq. (6)] has to be solved, which is however computationally less expensive than the other ones. 

As far as the one-dimensional test-case is concerned, the singlephase model of Wu and Wang (2007) and the two-phase model of Greco et al. (2012a) have been considered for comparison. The one-dimensional model by Wu and Wang (2007) is a single-phase mixture model, which considers both the suspended and bed load and accounts for variable bed load concentration. It is slightly less computationally expensive than the presented model, since it requires the solution of four differential equations, plus the bed evolution one. The inertia of the bed load sediment is considered through an empirical spatial lag between the actual bed load solid transport rate and the capacity value. As a consequence, in addition to the Manning coefficient, two empirical parameters defining the nonequilibrium adaptation length of total load sediment transport have to be defined. Moreover, a correction factor for the transport stage number in the van Rijn (1984) formula $( k _ { t } )$ is introduced. It has been shown by the authors that while the results’ sensitivity to the adaptation length value was limited, the correction factor $k _ { t }$ significantly affected the predicted erosion magnitude. The two-phase model of Greco et al. (2012a) is constituted by four conservation laws plus the bed deformation equation. The suspended sediment motion is not accounted for and the sediment concentration in the bed load is assumed to be constant. The concrete model application needs the estimation of the Chezy coefficient and of the bed load concentration. The latter has been assumed to be equal to the bed concentration (Greco et al. 2012a). 

Fig. 15 compares the results for the one-dimensional test of the proposed model and of the two considered literature models. Fig. 15 indicates an evident improvement of the present model with respect to that by Greco et al. (2012a). In particular, the latter model fails to reproduce the observed weak hydraulic jump, with a gradual variation of the free surface and a very different position of the downstream waterfront. A significant underestimation of the bed scour is also noted. Present results support the consideration formulated by Li et al. (2013) that the assumption of a constant bed load concentration may fail during highly unsteady flows. Conversely, the present model performs similarly to the mixture model by Wu and Wang (2007), both in terms of bottom elevation and free surface profile (Fig. 15). Although the mixture model may appear more attractive for the lower computational complexity, the agreement in the bed erosion significantly depends on the calibrated value of the correction factor k<sub>t</sub>. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/2c7b6ec4a914b2b602052972c75e9ff287c6863a8e4148e0118f537d40b98fd2.jpg)



(a)



x (m)



(b)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/47beefa053ee0999dac1732053f76c3a36a859caf6995823487c893bb979aad4.jpg)



x (m)



Fig. 15. One-dimensional dry-bed test; comparison with results from previous models: bottom and free surface profile: (a) $t = 0 . 4$ s; (b) $t = 0 . 5$ s after dam removal


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/296a58087cc2d946612d42fffeb5136d2eace4906714b590f18b3d15b3b74aa8.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/7a7195f1771d28922bd6545b03e27cfdd55f91350d19874e0ed7528b8184d4dd.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/3fc40f9309197becbb2fe9bb4e9e5eb89f4686d58258f8f8300e61d9a9febfa1.jpg)



(d)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/824a5381d1a2eeed3cd63990a9e633ff9b27c8cc96b29f022c103331c81fc2b4.jpg)



t (s)



Fig. 16. Two-dimensional dry-bed test; time series of free-surface elevation compared with results from previous models: (a) gauges 1 and 4; (b) gauges 2 and 3; (c) gauges 5 and 8; (d) gauges 6 and 7


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/99a2e3027144a0046015500670759352bce66bd74c149b32f5d913e141069f5c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/dfec5dd899b58dc4a1e7e2cb01896e2754b4e9565045a070b19c17de5efad790.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/70d56522479fc5d4f82998ebaeac620818cabe202cc81eb6e2432a1db65f6b64.jpg)



Fig. 17. Two-dimensional dry-bed test; final bottom profiles compared with results from previous models: (a) $y = 0 . 2 \ \mathrm { m } ;$ (b) $y = 0 . 7$ m; (c) y 1.45 m


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/4b8c529542f6ad3c825901fe649acfd48c451cf6fae6a75b9ce88a21c879e99f.jpg)



(a)



t (s)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/4103a42903b487425dbebf1866528ec6ac442b420065b06c87aa260eb5962930.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/011d676c805ce89b5ddf421d189e3f10394ce4681e5aa5d49da7fa743095f07b.jpg)



t (s)



(b)



t (s)



(c)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/bebf03b4bfba195a17ef417b8e4f841c21ddb3e282f19d43e0e5f7306da13630.jpg)



(d)



t (s)



Fig. 18. Two-dimensional wet-bed test; time series of free-surface elevation compared with results from previous models: (a) gauges 1 and 4; (b) gauges 2 and 3; (c) gauges 5 and 8; (d) gauges 6 and 7


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/05c3bfb0991da153126011738daccee45e675b0310749429d26ed878297b908f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/0234a3c471fa82e9a70be5d1b370f96c01c1d668899b2e22d7ce2ec0ab695604.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-24/c562ee27-55f2-48b7-8277-fbc4c1e288e6/8c7b3d646b9f0a3fa7b3e063346b0600a8fc74deb801e2143b336fb16f911eb8.jpg)



Fig. 19. Two-dimensional wet-bed test; final bottom profiles compared with results from previous models: (a) $y = 0 . 2$ m; (b) $y = 0 . 7$ m; (c) y 1.45 m


For the two-dimensional test cases, the comparison involves the single-phase model of Canelas et al. (2013) and the two-layer model of Swartenbroekx et al. (2013). The mixture two-dimensional model of Canelas et al. (2013) exhibits a much smaller computational complexity than the present one, being constituted by four conservation type laws plus the bed evolution one. Similar to the Wu and Wang (2007) model, a spatial lag between the actual bed load discharge and the equilibrium value is introduced to mimic the effects of the bed load inertia in the layer. The spatial lag is computed through an ad hoc formula that includes three additional calibration parameters fixed through a heuristic adjustment process. The computational complexity of the two-layer model of 

Swartenbroekx et al. (2013) is slightly smaller than that of the present model. Indeed, it is composed of six conservation equations plus the bed evolution one. Similarly to the two-phase model of Greco et al. (2012a, b), it does not account for the suspended load and the sediment concentration in the bed load is assumed constant. The sediment inertia in the bed load layer is fully described through the balance equation for the mixture momentum in the transport layer. The shear stresses between the layers are expressed through two constant friction factors, which have been determined through calibration against experimental results. 

Fig. 16 (resp. Fig. 18) compares the results of the present model for the two dimensional Test Case 1 (resp. Case 2) in terms of free-surface elevation with those of Canelas et al. (2013) and Swartenbroekx et al. (2013). Fig. 17 (resp. Fig. 19) is the counterpart of Fig. 16 (resp. Fig. 18) in terms of final topography. Both free-surface elevation history (Figs. 16 and 18) and final bottom topography (Figs. 17 and 19) are reproduced with an accuracy comparable to that of the model by Swartenbroekx et al. (2013) and with a slight improvement with respect to the mixture model of Canelas et al. (2013), despite the proper calibration of the three additional parameters. However, all models exhibit a slight but systematic under-prediction of the experimentally observed deposition. 

## Conclusions

A two-phase depth-averaged model able to deal with both bed load and suspended sediment transport has been proposed. The mathematical model, based on mass and momentum conservation equations for liquid and sediment phases, accounts for variable concentration both in the bed load and in the suspended load region. The entrainment/deposition of sediments from the bed toward the bed load layer is evaluated by a formula based on a modified van Rijn mobility parameter, whereas for the exchange between bed and suspended load a first-order exchange law is considered. The adopted set of closure relations is shown to comply, under uniform conditions of flow, with several empirical scaling laws for sediment transport and to allow for relatively accurate evaluation of the bed load layer thickness up to values of the Shields parameter order of unity. Two of the three dimensionless parameters of the model, the Chezy and the Meyer-Peter and Müller formula coefficients, may be evaluated based on extensive literature indications. The third one, $^ { c _ { 1 } , }$ is allowed to vary in a range limited by theoretically deduced lower and upper bounds. 

It has been proved that the proposed model is hyperbolic and the analytical expression of the eigenvalues has been provided. A numerical method based on a finite-volume approach has been used for the simulation of three experiments concerning three different dam breaks, showing a good agreement between simulated and experimental results. The results show that accounting for the variability concentration in the two-phase formulation leads to a neat improvement of the model performance. Finally, for all test, it has been demonstrated that the value of the free parameter $c _ { 1 }$ has only a marginal influence on the results’ quality. A further confirmation of this conclusion could be obtained through future application of the model to a wider class of morphodynamic transients. 

## References



Ancey, C., Andreini, N., and Epely-Chauvin, G. (2012). “Viscoplastic dambreak waves: Review of simple computational approaches and comparison with experiments.” Adv. Water Resour., 48, 79–91. 





Armanini, A. (2013). “Granular flows driven by gravity.” J. Hydraul. Res., 51(2), 111–120. 





Armanini, A., and Di Silvio, G. (1988). “A one-dimensional model for the transport of a sediment mixture in non-equilibrium conditions.” J. Hydraul. Res., 26(3), 275–292. 





Bagnold, R. A. (1956). “The flow of cohesionless grains in fluids.” Phil. Trans. Roy. Soc. A, 249(964), 235–297. 





Baldock, T. E., Tomkins, M. R., Nielsen, P., and Hughes, M. G. (2004). “Settling velocity of sediments at high concentrations.” Coast. Eng., 51(1), 91–100. 





Bridge, J. S., and Dominic, D. F. (1984). “Bed load grain velocities and sediment transport rates.” Water Resour. Res., 20(4), 476–490. 





Brooks, G. R., and Lawrence, D. E. (1999). “The drainage of the Lake Ha! Ha! reservoir and downstream impacts along Ha! Ha! River, Saguenay area, Quebec, Canada.” Geomorphol., 28(1–2), 141–167. 





Byrd, T. C., and Furbish, D. J. (2000). “Magnitude of deviatoric terms in vertically averaged flow equations.” Earth Surf. Processes Landforms, 25(3), 319–328. 





Canelas, R., Murillo, J., and Ferreira, R. M. (2013). “Two-dimensional depth-averaged modelling of dam-break flows over mobile beds.” J. Hydraul. Res., 51(4), 392–407. 





Cao, Z., and Carling, P. A. (2002). “Mathematical modelling of alluvial rivers: Reality and myth. Part 1: General overview.” Marit. Eng., 154(3), 207–219. 





Cao, Z., Day, R., and Egashira, S. (2002). “Coupled and decoupled numerical modeling of flow and morphological evolution in alluvial rivers.” J. Hydraul. Eng., 10.1061/(ASCE)0733-9429(2002)128:3(306), 306–321. 





Capart, H., and Young, D. (2002). “Two-layer shallow water computations of torrential geomorphic flows.” Proc., River Flow 2002, Swets & Zeitlinger, Lisse, Netherlands, 1003–1012. 





Capart, H., and Young, D. L. (1998). “Formation of a jump by the dambreak wave over a granular bed.” J. Fluid Mech., 372, 165–187. 





Chen, S. C., and Peng, S. H. (2006). “Two-dimensional numerical model of two-layer shallow water equations for confluence simulation.” Adv. Water Resour., 29(11), 1608–1617. 





Cheng, N.-S. (2006). “Influence of shear stress fluctuation on bed particle mobility.” Phys. Fluids, 18, 096602. 





Courant, R., and Hilbert, D. (1961). Methods of mathematical physics, Vol. 2, Wiley, New York. 





Defina, A., and Bixio, A. C. (2005). “Mean flow and turbulence in vegetated open channel flow.” Water Resour. Res., 41, W07006. 





Dewals, B., Rulot, F., Erpicum, S., Archambeau, P., and Pirotton, M. (2011). “Advanced topics in sediment transport modelling: Non-alluvial beds and hyperconcentrated flows, sediment transport.” 〈http://www .intechopen.com/books/sediment-transport/advanced-topics-in-sediment -transport-modelling-non-alluvial-beds-and-hyperconcentrated-flows〉 (Jan. 12, 2014). 





Di Cristo, C., Evangelista, S., Leopardi, A., Greco, M., and Iervolino, M. (2014a). “Numerical simulation of a dam-break with a wide range of shields parameter.” Proc., Int. Conf. on Fluvial Hydraulics, River Flow 2014, CRC, 1680–1687. 





Di Cristo, C., Iervolino, M., and Vacca, A. (2006). “Linear stability analysis of a 1-D model with dynamical description of a bed load transport.” J. Hydraul. Res., 44(4), 480–487. 





Di Cristo, C., Iervolino, M., and Vacca, A. (2014b). “Applicability of kinematic, diffusion and quasi-steady dynamic wave models to shallow mud flows.” J. Hydrol. Eng., 10.1061/(ASCE)HE.1943-5584.0000881, 956–965. 





Di Cristo, C., Iervolino, M., and Vacca, A. (2014c). “Simplified wave models applicability to shallow mud flows modeled as power-law fluids.” J. Mt. Sci., 11(6), 1454–1465. 





Di Cristo, C., Iervolino, M., and Vacca, A. (2015). “Diffusive approximation for unsteady mud flows with backwater effect.” Adv. Water Resour., 81, 84–94. 





Duran, O., Andreotti, B., and Claudin, P. (2012). “Numerical simulation of turbulent sediment transport, from bed load to saltation.” Phys. Fluids, 24(10), 1–23. 





El Kadi Abderrezzak, K., and Paquier, A. (2011). “Applicability of sediment transport capacity formulas to dam-break flows over movable beds.” J. Hydraul. Eng., 10.1061/(ASCE)HY.1943-7900.0000298, 209–221. 





Evangelista, S., Altinakar, M. S., Di Cristo, C., and Leopardi, A. (2013). “Simulation of dam-break waves on movable beds using a multi-stage centred scheme.” Int. J. Sediment. Res., 28(3), 269–284. 





Fernandez-Luque, R., and Van Beek, R. (1976). “Erosion and transport of bed-load sediment.” J. Hydraul. Res., 14(2), 127–144. 





Fraccarollo, L., and Capart, H. (2002). “Riemann wave description of erosional dam-break flows.” J. Fluid Mech., 461, 183–228. 





Furbish, D. J., Haff, P. K., Roseberry, J. C., and Schmeeckle, M. W. (2012). “A probabilistic description of the bed load sediment flux: 1. Theory.” J. Geophys. Res., 117(F3), F03031. 





Garegnani, G., Rosatti, G., and Bonaventura, L. (2011). “Free surface flows over mobile bed: Mathematical analysis and numerical modeling of coupled and decoupled approaches.” Commun. Appl. Ind. Math., 2(1), 1–22. 





Graf, W. H. (1998). Fluvial hydraulics: Flow and transport processes in channels of simple geometry, Wiley, Chichester, U.K. 





Greco, M., Iervolino, M., and Leopardi, A. (2008a). “Discussion on divergence form for bed slope source term in shallow water equations.” J. Hydraul. Eng., 10.1061/(ASCE)0733-9429(2008)134:5(676), 676–678. 





Greco, M., Iervolino, M., Leopardi, A., and Vacca, A. (2012a). “A twophase model for fast geomorphic shallow flows.” Int. J. Sediment. Res., 27(4), 409–425. 





Greco, M., Iervolino, M., and Vacca, A. (2008b). “Discussion on boundary conditions in a two-layer geomorphological model: Application to a hydraulic jump over a mobile bed.” J. Hydraul. Res., 46(6), 856–860. 





Greco, M., Iervolino, M., Vacca, A., and Leopardi, A. (2012b). “Two-phase modelling of total sediment load in fast geomorphic transients.” River Flow 2012, Proc., Int. Conf. on Fluvial Hydraulics, Vol. 1, Colegio de Ingenieros Civiles de Costa Rica (CiC), 643–648. 





Griffiths, G. A. (1981). “Flow resistance in coarse gravel bed rivers.” J. Hydraul. Div., 107(7), 899–918. 





Julien, P. Y., and Bounvilay, B. (2013). “Velocity of rolling bed load particles.” J. Hydraul. Eng., 10.1061/(ASCE)HY.1943-7900.0000657, 177–186. 





Keylock, C. J., Hardy, R. J., Parsons, D. R., Ferguson, R. I., Lane, S. N., and Richards, K. S. (2005). “The theoretical foundations and potential for large-eddy simulation (LES) in fluvial geomorphic and sedimentological research.” Earth-Sci. Rev., 71(3–4), 271–304. 





Lajeunesse, E., Malverti, L., and Charru, F. (2010). “Bed load transport in turbulent flow at the grain scale: Experiments and modelling.” J. Geophys. Res. F Earth Surface, 115(F4), F04001. 





Lamb, M. P., Dietrich, W. E., and Venditti, J. G. (2008). “Is the critical Shields stress for incipient sediment motion dependent on channelbed slope?” J. Geophys. Res., 113(F2), 1–23. 





Lee, H. Y., Chen, Y. H., You, J. Y., and Lin, Y. T. (2000). “Investigations of continuous bed load saltating process.” J. Hydraul. Eng., 10.1061/ (ASCE)0733-9429(2000)126:9(691), 691–700. 





Lee, H. Y., and Hsu, I. S. (1994). “Investigation of saltating particle motion.” J. Hydraul. Eng., 10.1061/(ASCE)0733-9429(1994)120:7(831), 831–845. 





Leopardi, A., Oliveri, E., and Greco, M. (2002). “Two-dimensional modeling of floods to map risk prone areas.” J. Water Resour. Plann. Manage., 10.1061/(ASCE)0733-9496(2002)128:3(168), 168–178. 





Li, J., Cao, Z., Pender, G., and Liu, Q. (2013). “A double layer-averaged model for dam-break flows over mobile bed.” J. Hydraul. Res., 51(5), 518–534. 





Lopez, F., and Garcia, M. H. (1996). “Turbulence structure in cobblebed open-channel flow.” Civil engineering studies, Univ. of Illinois, Urbana, IL. 





Marsooli, R., and Wu, W. (2015). “Three-dimensional numerical modeling of dam-break flows with sediment transport over movable beds.” J. Hydraul. Eng.,10.1061/(ASCE)HY.1943-7900.0000947, 04014066. 





Maude, A. D., and Whitmore, R. L. (1958). “A generalized theory of sedimentation.” Br. J. Appl., 9(12), 477–482. 





Meland, N., and Norrman, J. O. (1966). “Transport velocities of single particles in bed load motion.” Geografiska Annaler, 48A(4), 165–182. 





Meyer-Peter, E., and Müller, R. (1948). “Formulas for bed-load transport.” Proc., Int. Association of Hydraulic Research 2nd Meeting, Stockholm. 





Nikora, V., and Goring, D. (2000). “Flow turbulence over fixed and weakly mobile gravel beds.” J. Hydraul. Eng., 10.1061/(ASCE)0733-9429 (2000)126:9(679), 679–690. 





Nino, Y., and Garcia, M. (1998). “Experiments on saltation of sand in water.” J. Hydraul. Eng., 10.1061/(ASCE)0733-9429(1998)124:10(1014), 1014–1025. 





Nino, Y., Garcia, M., and Ayala, L. (1994). “Gravel saltation. 1. Experiments.” Water Resour. Res., 30(6), 1907–1914. 





Papanicolaou, A. N., Diplas, P., Evaggelopoulos, N., and Fotopoulos, S. (2002a). “Stochastic incipient motion criterion for spheres under various bed packing conditions.” J. Hydraul. Eng., 10.1061/(ASCE)0733 -9429(2002)128:4(369), 369–380. 





Papanicolaou, A. N., Knapp, D., and Strom, K. (2002b). “Bedload predictions by using the concept of particle velocity: Applications.” Proc., ASCE/EWRI and IAHR Int. Conf. on Hydraulic Measurements and Experimental Methods, ASCE, Reston, VA, 1–10. 





Parker, G., Seminara, G., and Solari, L. (2003). “Bed load at low Shields stress on arbitrarily sloping beds: Alternative entrainment formulation.” Water Resour. Res., 39(7), 1183–1194. 





Pelanti, M., Bouchut, F., and Mangeney, A. (2008). “A Roe-type scheme for two-phase shallow granular flows over variable topography.” Math. Model. Numer. Anal., 42(5), 851–885. 





Pitman, E. B., and Le, L. (2005). “A two-fluid model for avalanche and debris flows.” Phil. Trans. R. Soc. A, 363(1832), 1573–1601. 





Pontillo, M., Schmocker, L., Greco, M., and Hager, W. H. (2010). “1D numerical evaluation of dike erosion due to overtopping.” J. Hydraul. Res., 48(5), 573–582. 





Pudasaini, S. P. (2012). “A general two-phase debris flow model.” J. Geophys. Res., 117(F3), F03010. 





Pudasaini, S. P., Wang, Y., and Hutter, K. (2005). “Modelling debris flows down general channels.” Nat. Hazards Earth Syst. Sci., 5(6), 799–819. 





Ramesh, B., Kothyari, U. C., and Murugesan, K. (2011). “Near-bed particle motion over transitionally-rough bed.” J. Hydraul. Res., 49(6), 757–765. 





Richardson, J. F., and Zaki, W. N. (1954). “Sedimentation and fluidisation: Part 1.” Trans. Inst. Chem. Eng., 32, 35–53. 





Rosatti, G., and Begnudelli, L. (2013). “A closure-independent generalized Roe solver for free-surface, two-phase flows over mobile bed.” J. Comput. Phys., 255, 362–383. 





Sabbagh-Yazdi, S., and Jamshidi, M. (2013). “Depth-averaged hydrodynamic model for gradual breaching of embankment dams attributable to overtopping considering suspended sediment transport.” J. Hydraul. Eng., 10.1061/(ASCE)HY.1943-7900.0000706, 580–592. 





Savary, C., and Zech, Y. (2007). “Boundary conditions in a two-layer geomorphological model: Application to a hydraulic jump over a mobile bed.” J. Hydraul. Res., 45(3), 316–332. 





Savary, C., and Zech, Y. (2008). “Boundary conditions in a two-layer geomorphological model: Application to a hydraulic jump over a mobile bed.” J. Hydraul. Res., 46(6), 858–860. 





Sekine, M., and Kikkawa, H. (1992). “Mechanics of saltating grains. II.” J. Hydraul. Eng., 10.1061/(ASCE)0733-9429(1992)118:4(513), 513–535. 





Seminara, G., Solari, L., and Parker, G. (2002). “Bed load at low Shields stress on arbitrarily sloping beds: Failure of the Bagnold hypothesis.” Water Resour. Res., 38(11), 1249–1271. 





Simpson, G., and Castelltort, S. (2006). “Coupled model of surface water flow, sediment transport and morphological evolution.” Comput. Geosci., 32(10), 1600–1614. 





Singh, A. K., Kothyari, U. C., and Ranga Raju, K. G. (2004). “Rapidly varying transient flows in alluvial rivers.” J. Hydraul. Res., 42(5), 473–486. 





Soares-Frazão, S. (2012). “Dam-break flows over mobile beds: Experiments and benchmark tests for numerical models.” J. Hydraul. Res., 50(4), 364–375. 





Soldati, A., and Marchioli, C. (2012). “Sediment transport in steady turbulent boundary layers: Potentials, limitations, and perspectives for Lagrangian tracking in DNS and LES.” Adv. Water Resour., 48, 18–30. 





Spinewine, B., and Zech, Y. (2007). “Small-scale laboratory dam-break waves on movable beds.” J. Hydraul. Res., 45(extra issue), 73–86. 





Sturm, T. W. (2013). “Hydraulic engineering: A rising wave of progress.” J. Hydraul. Eng., 10.1061/(ASCE)HY.1943-7900.0000689, 111–113. 





Swartenbroekx, C., Zech, Y., and Soares-Frazão, S. (2013). “Twodimensional two-layer shallow water model for dam break flows with significant bed load transport.” Int. J. Numer. Meth. Fluids, 73(5), 477–508. 





Valiani, A., and Begnudelli, L., (2006). “Divergence form for bed slope source term in shallow water equations.” J. Hydraul. Eng., 10.1061/ (ASCE)0733-9429(2006)132:7(652), 652–665. 





van Rijn, L. C. (1984). “Sediment pick-up functions.” J. Hydraul. Eng., 10.1061/(ASCE)0733-9429(1984)110:10(1494), 1494–1502. 





Wang, S. S. Y., and Wu, W. M. (2005). “Computational simulation of river sedimentation and morphology—A review of the state of the art.” Int. J. Sediment Res., 20(1), 7–29. 





Wiberg, P. L., and Smith, J. D. (1985). “A theoretical model for saltating grains in water.” J. Geophys. Res., 90(c4), 7341–7354. 





Wohl, E. E., and Thompson, D. M. (2000). “Velocity characteristics along a small step-pool channel.” Earth Surf. Processes Landforms, 25(4), 353–367. 





Wong, M., and Parker, G. (2006). “Reanalysis and correction of bed-load relation of Meyer-Peter and Müller using their own database.” J. Hydraul. Eng., 10.1061/(ASCE)0733-9429(2006)132:11(1159), 1159–1168. 





Wong, M., Parker, G., De Vries, P., Brown, T. M., and Burges, S. J. (2007). “Experiments on dispersion of tracer stones under lower-regime planebed equilibrium bed load transport.” Water Resour. Res., 43(3), 1–23. 





Wu, F. C., and Chou, Y. J. (2003). “Rolling and lifting probabilities for sediment entrainment.” J. Hydraul. Eng., 10.1061/(ASCE)0733-9429 (2003)129:2(110), 110–119. 





Wu, W. (2007). Computational river dynamics, Taylor & Francis, London. 





Wu, W., and Wang, S. S.-Y. (2007). “One dimensional modeling of dambreak flow over movable beds.” J. Hydraul. Eng., 10.1061/(ASCE) 0733-9429(2007)133:1(48), 48–58. 





Wu, W., Wang, S. S.-Y., and Jia, Y. (2000). “Nonuniform sediment transport in alluvial rivers.” J. Hydraul. Res., 38(6), 427–434. 





Xia, J., Lin, B., Falconer, R. A., and Wang, G. (2010). “Modelling dambreak flows over mobile beds using a 2D coupled approach.” Adv. Water Resour., 33(2), 171–183. 





Zhang, S., Duan, J. G., and Strelkoff, T. S. (2013). “Grain-scale nonequilibrium sediment-transport model for unsteady flow.” J. Hydraul. Eng., 10.1061/(ASCE)HY.1943-7900.0000645, 22–36. 

