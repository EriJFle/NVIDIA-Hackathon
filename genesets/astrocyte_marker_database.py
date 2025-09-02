# Core astrocyte markers (literature-validated panel)
# Based on Batiuk et al. 2020, Bayraktar et al. 2020, and Zeisel et al. 2018
astrocyte_markers = [
    # Pan-astrocyte core markers
    "Aldh1l1",   # Aldehyde dehydrogenase 1 family, member L1 - gold standard
    "Slc1a3",    # Solute carrier family 1, member 3 (GLAST) - glutamate transporter
    "Aqp4",      # Aquaporin 4 - water channel, blood-brain barrier
    "Gja1",      # Gap junction protein, alpha 1 (Connexin 43)
    "Glul",      # Glutamate-ammonia ligase (glutamine synthetase)
    
    # Additional pan-astrocyte markers
    "Sox9",      # SRY-box 9 - transcription factor, astrocyte identity
    "Nfia",      # Nuclear factor I/A - astrocyte development
    "Fgfr3",     # Fibroblast growth factor receptor 3
    "Kcnj10",    # Potassium channel, inwardly rectifying subfamily J, member 10 (Kir4.1)
    "Gjb6",      # Gap junction protein, beta 6 (Connexin 30)
    "Slc1a2",    # Solute carrier family 1, member 2 (GLT-1) - main glutamate transporter
    "S100b",     # S100 calcium binding protein B
    "Gfap",      # Glial fibrillary acidic protein - classical but not universal
    "Vim",       # Vimentin - intermediate filament
    "Fabp7",     # Fatty acid binding protein 7 (brain lipid binding protein)
    "Slc4a4",    # Solute carrier family 4, member 4 (NBC1) - bicarbonate transporter
    "Mlc1",      # Megalencephalic leukoencephalopathy with subcortical cysts 1
    "Cldn10",    # Claudin 10 - tight junction
]

# Pan-reactive astrocyte markers (A1/A2 and general reactivity)
# Based on Liddelow et al. 2017, Zamanian et al. 2012, and recent single-cell studies
pan_reactive_astrocyte_markers = [
    # A1-like (neuroinflammatory) reactive markers
    "Lcn2",      # Lipocalin 2 - iron sequestration, inflammation
    "Serpina3n", # Serpin peptidase inhibitor, clade A, member 3N
    "C3",        # Complement component 3
    "H2-D1",     # Histocompatibility 2, D region locus 1 (MHC class I)
    "Gbp2",      # Guanylate binding protein 2
    "Psmb8",     # Proteasome subunit, beta type, 8
    "Iigp1",     # Interferon inducible GTPase 1
    "Srgn",      # Serglycin
    "Amigo2",    # Adhesion molecule with Ig-like domain 2
    
    # A2-like (neuroprotective) reactive markers
    "S100a10",   # S100 calcium binding protein A10
    "Cd109",     # CD109 molecule
    "Emp1",      # Epithelial membrane protein 1
    "Tm4sf1",    # Transmembrane 4 L six family member 1
    "B3gnt5",    # UDP-GlcNAc:betaGal beta-1,3-N-acetylglucosaminyltransferase 5
    "Clcf1",     # Cardiotrophin-like cytokine factor 1
    "Ptx3",      # Pentraxin 3
    "Tgm1",      # Transglutaminase 1
    "Sphk1",     # Sphingosine kinase 1
    
    # General reactive markers
    "Gfap",      # Glial fibrillary acidic protein - upregulated in reactivity
    "Vim",       # Vimentin - upregulated in reactivity
    "Timp1",     # TIMP metallopeptidase inhibitor 1
    "Cxcl10",    # Chemokine (C-X-C motif) ligand 10
    "Osmr",      # Oncostatin M receptor
    "Cd44",      # CD44 molecule - hyaluronic acid receptor
    "Hspb1",     # Heat shock protein family, member 1
    "Aspg",      # Asparaginase
    "Steap4",    # STEAP family member 4
    "Cp",        # Ceruloplasmin
]

# Region-specific astrocyte markers
# Based on Batiuk et al. 2020 and Bayraktar et al. 2020
cortical_astrocyte_markers = [
    "Gpc5",      # Glypican 5 - cortical layers
    "Luzp2",     # Leucine zipper protein 2
    "Dner",      # Delta/notch-like EGF repeat containing
    "Ntsr2",     # Neurotensin receptor 2
    "Grm5",      # Glutamate receptor, metabotropic 5
]

hippocampal_astrocyte_markers = [
    "Gfap",      # Higher in hippocampal astrocytes
    "Aqp4",      # Enriched in hippocampus
    "Aldoc",     # Aldolase C, fructose-bisphosphate
    "Prox1",     # Prospero homeobox 1
]

white_matter_astrocyte_markers = [
    "Gfap",      # Fibrous astrocytes
    "Cspg4",     # Chondroitin sulfate proteoglycan 4
    "Pdgfra",    # Platelet derived growth factor receptor, alpha polypeptide
]

# Developmental astrocyte markers
developmental_astrocyte_markers = [
    "Sox9",      # Early astrocyte specification
    "Nfia",      # Astrocyte differentiation
    "Nfib",      # Nuclear factor I/B
    "Hes5",      # Hes family bHLH transcription factor 5
    "Id3",       # Inhibitor of DNA binding 3
    "Fabp7",     # Radial glia/early astrocytes
    "Nestin",    # Nestin - neural stem cell marker
]

# Disease-associated astrocyte markers
disease_astrocyte_markers = [
    # Alzheimer's disease
    "Gfap",      # Upregulated in AD
    "S100b",     # Upregulated in AD
    "C3",        # Complement activation
    "Cxcl10",    # Chemokine signaling
    
    # Multiple sclerosis/demyelination
    "Lcn2",      # Iron handling
    "Cxcl10",    # T cell recruitment
    "Ccl2",      # Chemokine (C-C motif) ligand 2
    "Il6",       # Interleukin 6
    
    # Epilepsy
    "Gfap",      # Reactive gliosis
    "Kcnj10",    # Potassium buffering impairment
    "Aqp4",      # Water homeostasis
    "S100b",     # Calcium signaling
    
    # Stroke/ischemia
    "Lcn2",      # Acute phase response
    "Hif1a",     # Hypoxia-inducible factor 1, alpha subunit
    "Vegfa",     # Vascular endothelial growth factor A
    "Hmox1",     # Heme oxygenase 1
]

# Excitatory neuron markers
excitatory_neuron_markers = [
    "Slc17a7", "Slc17a6", "Tbr1", "Satb2",
    "Cux1", "Cux2", "Rorb", "Bcl11b",
]

# Inhibitory neuron markers
inhibitory_neuron_markers = [
    "Gad1", "Gad2", "Slc32a1", "Pvalb",
    "Sst", "Vip", "Reln",
]

# OPC markers
opc_markers = [
    "Pdgfra", "Cspg4", "Sox10", "Olig1",
    "Olig2", "Gpr17", "Ptprz1", "Vcan",
]

# Oligodendrocyte markers
oligodendrocyte_markers = [
    "Mbp", "Plp1", "Mog", "Cldn11",
    "Mag", "Mobp", "Opalin",
]

# Microglia markers
microglia_markers = [
    "Tmem119", "P2ry12", "Hexb", "Cx3cr1",
    "Trem2", "Tyrobp",
]

# Endothelial cell markers
endothelial_markers = [
    "Pecam1", "Cdh5", "Vwf", "Cldn5",
    "Ocln", "Kdr", "Esam", "Erg",
]


pericyte_markers = [
    "Pdgfrb", "Rgs5", "Abcc9", "Kcnj8",
    "Notch3", "Higd1b", "Ndufa4l2", "Kcnj13",
]

# Vascular smooth muscle cell markers
smooth_muscle_cell_markers = [
    "Acta2", "Myh11", "Tagln", "Cnn1",
    "Myl9", "Lmod1",
]

# Ependymal cell markers
ependymal_markers = [
    "Foxj1", "Rsph1", "Pifo", "Dnah5",
    "Tekt1", "Dynlrb2",
]

# Tanycyte markers
tanycyte_markers = [
    "Rax", "Lhx2", "Dio2", "Gpr50",
    "Slc16a2", "Ptn",
]

# Choroid plexus markers
choroid_plexus_markers = [
    "Ttr", "Folr1", "Aqp1", "Otx2",
    "Kcne2",
]

# Comprehensive marker dictionary
comprehensive_markers = {
    # Astrocyte markers
    "astrocyte_core": astrocyte_markers,
    "astrocyte_reactive": pan_reactive_astrocyte_markers,
    "astrocyte_cortical": cortical_astrocyte_markers,
    "astrocyte_hippocampal": hippocampal_astrocyte_markers,
    "astrocyte_white_matter": white_matter_astrocyte_markers,
    "astrocyte_developmental": developmental_astrocyte_markers,
    "astrocyte_disease": disease_astrocyte_markers,
    
    # Non-astrocyte markers
    "excitatory_neuron": excitatory_neuron_markers,
    "inhibitory_neuron": inhibitory_neuron_markers,
    "opc": opc_markers,
    "oligodendrocyte": oligodendrocyte_markers,
    "microglia": microglia_markers,
    "endothelial": endothelial_markers,
    "pericyte": pericyte_markers,
    "smooth_muscle_cell": smooth_muscle_cell_markers,
    "ependymal": ependymal_markers,
    "tanycyte": tanycyte_markers,
    "choroid_plexus": choroid_plexus_markers,
}

# Backward compatibility - minimal markers
minimal_markers = {
    "astrocyte": astrocyte_markers,
    "pan_reactive_astrocyte": pan_reactive_astrocyte_markers,
    "excitatory_neuron": excitatory_neuron_markers,
    "inhibitory_neuron": inhibitory_neuron_markers,
    "opc": opc_markers,
    "oligodendrocyte": oligodendrocyte_markers,
    "microglia": microglia_markers,
    "endothelial": endothelial_markers,
    "pericyte": pericyte_markers,
    "smooth_muscle_cell": smooth_muscle_cell_markers,
    "ependymal": ependymal_markers,
    "tanycyte": tanycyte_markers,
    "choroid_plexus": choroid_plexus_markers,
}

# Exclusion markers (all non-astrocyte cell types)
exclude_markers = (excitatory_neuron_markers + inhibitory_neuron_markers + 
                  opc_markers + oligodendrocyte_markers + microglia_markers + 
                  endothelial_markers + pericyte_markers + smooth_muscle_cell_markers + 
                  ependymal_markers + tanycyte_markers + choroid_plexus_markers)

# Functions for marker retrieval
def get_astrocyte_markers(region=None, state=None, development=None):
    """
    Get astrocyte markers based on specific criteria.
    
    Parameters
    ----------
    region : str, optional
        Specific brain region ('cortical', 'hippocampal', 'white_matter')
    state : str, optional  
        Astrocyte state ('reactive', 'disease')
    development : bool, optional
        Include developmental markers
        
    Returns
    -------
    list
        List of relevant astrocyte markers
    """
    markers = astrocyte_markers.copy()
    
    if region == 'cortical':
        markers.extend(cortical_astrocyte_markers)
    elif region == 'hippocampal':
        markers.extend(hippocampal_astrocyte_markers)
    elif region == 'white_matter':
        markers.extend(white_matter_astrocyte_markers)
    
    if state == 'reactive':
        markers.extend(pan_reactive_astrocyte_markers)
    elif state == 'disease':
        markers.extend(disease_astrocyte_markers)
    
    if development:
        markers.extend(developmental_astrocyte_markers)
    
    return list(set(markers))  # Remove duplicates

def get_contamination_markers():
    """Get all non-astrocyte markers for contamination detection."""
    return exclude_markers

def get_marker_categories():
    """Get all available marker categories."""
    return list(comprehensive_markers.keys())