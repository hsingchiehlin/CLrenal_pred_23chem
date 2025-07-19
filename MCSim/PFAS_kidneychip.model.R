# -------------------------------------------------------------------------
# PFAS_kidneychip.model
#
# Units:
# - time in hours
# - volumes in liters
# - masses of substances in micromoles
# - concentrations of substances in microM
# -------------------------------------------------------------------------

States  = { 
  
  # 1: Lumen (A) to Blood (B) study
  A1_blood,      # Quantity in blood compartment (micromoles)
  A1_cell,       # Quantity in cell lysate compartment (micromoles)
  A1_lumen,      # Quantity in lumen compartment (micromoles)
  
  # 2: Blood (B) to Lumen (A) study
  A2_blood,      # Quantity in blood compartment (micromoles)
  A2_cell,       # Quantity in cell lysate compartment (micromoles)
  A2_lumen,      # Quantity in lumen compartment (micromoles)
  
  # 3: Both add study
  #A3_blood,      # Quantity in blood compartment (micromoles)
  #A3_cell,       # Quantity in cell lysate compartment (micromoles)
  #A3_lumen       # Quantity in lumen compartment (micromoles)
  
};


Outputs = { 

  # 1: Lumen (A) to Blood (B) study
  C1_blood,      # concentration in blood compartment (uM)
  C1_cell,       # concentration in cell lysate compartment (uM)
  C1_lumen,      # concentration in lumen compartment (uM)
  
  # 2: Blood (B) to Lumen (A) study
  C2_blood,      # concentration in blood compartment (uM)
  C2_cell,       # concentration in cell lysate compartment (uM)
  C2_lumen,      # concentration in lumen compartment (uM)
  
  # 3: Both add study
  #C3_blood,      # Quantity in blood compartment (uM)
  #C3_cell,       # Quantity in cell lysate compartment (uM)
  #C3_lumen       # Quantity in lumen compartment (uM)
  
};         

Inputs = { 
  
  A1_lumen_init,
  A2_blood_init,
  #A3_blood_init,
  #A3_lumen_init
  
};       

p_1;             # cm/hr
p_2;             # cm/hr
p_3;             # cm/hr
p_4;             # cm/hr
area = 0.33;                # cm^2
V_blood = 0.00020;          # 0.20*0.001; 0.2ml -> L
V_cell = 5.748817e-8;       # (0.00065)^3*3.14*(4/3)*50000/1000
V_lumen = 0.00075;          # 0.75ml*0.001 -> L
GSD_error_C1_blood;
GSD_error_C1_cell;
GSD_error_C1_lumen;
GSD_error_C2_blood;
GSD_error_C2_cell;
GSD_error_C2_lumen;
GSD_error_C3_cell;
P1;
P2;
P3;
P4;

Initialize { 
  
  P1 = exp(p_1);
  P2 = exp(p_2);
  P3 = exp(p_3);
  P4 = exp(p_4);

};


Dynamics { 
  
  C1_blood = A1_blood/V_blood;
  C1_cell = A1_cell/V_cell;
  C1_lumen = A1_lumen/V_lumen;
  
  C2_blood = A2_blood/V_blood;
  C2_cell = A2_cell/V_cell;
  C2_lumen = A2_lumen/V_lumen;
  
  #C3_blood = A3_blood/V_blood;
  #C3_cell = A3_cell/V_cell;
  #C3_lumen = A3_lumen/V_lumen;
  
  #C1_blood   = (A1_blood  > 0 ? A1_blood/V_blood  : 1e-30); 
  #C1_cell    = (A1_cell   > 0 ? A1_cell/V_cell    : 1e-30); 
  #C1_lumen   = (A1_lumen  > 0 ? A1_lumen/V_lumen  : 1e-30); 
  
  #C2_blood   = (A2_blood  > 0 ? A2_blood/V_blood  : 1e-30); 
  #C2_cell    = (A2_cell   > 0 ? A2_cell/V_cell    : 1e-30); 
  #C2_lumen   = (A2_lumen  > 0 ? A2_lumen/V_lumen  : 1e-30); 
  
  #C3_blood   = (A3_blood  > 0 ? A3_blood/V_blood  : 1e-30); 
  #C3_cell    = (A3_cell   > 0 ? A3_cell/V_cell    : 1e-30); 
  #C3_lumen   = (A3_lumen  > 0 ? A3_lumen/V_lumen  : 1e-30); 
  
  
  #dt (A1_blood) = (A1_blood  > 0 ? (P1*C1_cell - P2*C1_blood)*area/1000  : 1e-30); 
  #dt (A1_cell)  = (A1_cell   > 0 ? (P2*C1_blood + P3*C1_lumen - (P1+P4)*C1_cell)*area/1000  : 1e-30);
  #dt (A1_lumen) = (A1_lumen  > 0 ? (P4*C1_cell - P3*C1_lumen)*area/1000  : 1e-30);
  
  #dt (A2_blood) = (A2_blood  > 0 ? (P1*C2_cell - P2*C2_blood)*area/1000  : 1e-30);
  #dt (A2_cell)  = (A2_cell   > 0 ? (P2*C2_blood + P3*C2_lumen - (P1+P4)*C2_cell)*area/1000  : 1e-30);
  #dt (A2_lumen) = (A2_lumen  > 0 ? (P4*C2_cell - P3*C2_lumen)*area/1000  : 1e-30);
  
  #dt (A3_blood) = (A3_blood  > 0 ? (P1*C3_cell - P2*C3_blood)*area/1000  : 1e-30);
  #dt (A3_cell)  = (A3_cell   > 0 ? (P2*C3_blood + P3*C3_lumen - (P1+P4)*C3_cell)*area/1000  : 1e-30);
  #dt (A3_lumen) = (A3_lumen  > 0 ? (P4*C3_cell - P3*C3_lumen)*area/1000  : 1e-30);
  
  dt (A1_blood) = (P1*C1_cell - P2*C1_blood)*area/1000;
  dt (A1_cell)  = (P2*C1_blood + P3*C1_lumen - (P1+P4)*C1_cell)*area/1000;
  dt (A1_lumen) = (P4*C1_cell - P3*C1_lumen)*area/1000;
  
  dt (A2_blood) = (P1*C2_cell - P2*C2_blood)*area/1000;
  dt (A2_cell)  = (P2*C2_blood + P3*C2_lumen - (P1+P4)*C2_cell)*area/1000;
  dt (A2_lumen) = (P4*C2_cell - P3*C2_lumen)*area/1000;
  
  #dt (A3_blood) = (P1*C3_cell - P2*C3_blood)*area/1000;
  #dt (A3_cell)  = (P2*C3_blood + P3*C3_lumen - (P1+P4)*C3_cell)*area/1000;
  #dt (A3_lumen) = (P4*C3_cell - P3*C3_lumen)*area/1000;
  
  #dt (A1_blood) = ((P1*C1_cell - P2*C1_blood)  > 0 ? (P1*C1_cell - P2*C1_blood)*area/1000  : 0); 
  #dt (A1_cell)  = ((P2*C1_blood + P3*C1_lumen - (P1+P4)*C1_cell)  > 0 ? (P2*C1_blood + P3*C1_lumen - (P1+P4)*C1_cell)*area/1000  : 0);
  #dt (A1_lumen) = ((P4*C1_cell - P3*C1_lumen)  > 0 ? (P4*C1_cell - P3*C1_lumen)*area/1000  : 0);
  
  #dt (A2_blood) = ((P1*C2_cell - P2*C2_blood)  > 0 ? (P1*C2_cell - P2*C2_blood)*area/1000  : 0);
  #dt (A2_cell)  = ((P2*C2_blood + P3*C2_lumen - (P1+P4)*C2_cell)   > 0 ? (P2*C2_blood + P3*C2_lumen - (P1+P4)*C2_cell)*area/1000  : 0);
  #dt (A2_lumen) = ((P4*C2_cell - P3*C2_lumen)  > 0 ? (P4*C2_cell - P3*C2_lumen)*area/1000  : 0);
  
  #dt (A3_blood) = ((P1*C3_cell - P2*C3_blood)  > 0 ? (P1*C3_cell - P2*C3_blood)*area/1000  : 0);
  #dt (A3_cell)  = ((P2*C3_blood + P3*C3_lumen - (P1+P4)*C3_cell)   > 0 ? (P2*C3_blood + P3*C3_lumen - (P1+P4)*C3_cell)*area/1000  : 0);
  #dt (A3_lumen) = ((P4*C3_cell - P3*C3_lumen)  > 0 ? (P4*C3_cell - P3*C3_lumen)*area/1000  : 0);
 
};

CalcOutputs { 
  
  C1_blood   = (C1_blood  > 0 ? A1_blood/V_blood  : 1e-15); 
  C1_cell    = (C1_cell   > 0 ? A1_cell/V_cell    : 1e-15); 
  C1_lumen   = (C1_lumen  > 0 ? A1_lumen/V_lumen  : 1e-15); 
  
  C2_blood   = (C2_blood  > 0 ? A2_blood/V_blood  : 1e-15); 
  C2_cell    = (C2_cell   > 0 ? A2_cell/V_cell    : 1e-15); 
  C2_lumen   = (C2_lumen  > 0 ? A2_lumen/V_lumen  : 1e-15); 
  
  
};

End.

