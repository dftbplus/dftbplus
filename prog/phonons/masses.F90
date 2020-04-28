module dftbp_masses
  use dftbp_accuracy
  implicit none
  private

  public getAtomicMass 

  contains 

  !!* Returns the atomic mass in a.m.u.
  !!* specie is a character string 
  function getAtomicMass(specie) result(mass)
    character(*), intent(in) :: specie
    real(dp) :: mass

    select case(trim(specie))
    case("H")
      mass = 1.00794_dp
    case("He")
      mass = 4.002602_dp
    case("Li")
      mass = 6.941_dp
    case("Be")
      mass = 9.012182_dp        
    case("B")
      mass = 10.811_dp 
    case("C")
      mass = 12.0107_dp
    case("N")
      mass = 14.00674_dp
    case("O")
      mass = 15.9994_dp
    case("F")
      mass = 18.9984032_dp
    case("Ne")
      mass = 20.1797_dp
    case("Na")
      mass = 22.989770_dp
    case("Mg")
      mass = 24.3050_dp
    case("Al")
      mass = 26.981538_dp
    case("Si")
      mass = 28.0855_dp
    case("P")
      mass = 30.973761_dp
    case("S")
      mass = 32.066_dp
    case("Cl")
      mass = 35.4527_dp
    case("Ar")
      mass = 39.948_dp
    case("K")
      mass = 39.0983_dp
    case("Ca")
      mass = 40.078_dp
    case("Sc")
      mass = 44.955910_dp
    case("Ti")
      mass = 47.867_dp     
    case("V")
      mass = 50.9415_dp    
    case("Cr")
      mass = 51.9961_dp    
    case("Mn")
      mass = 54.938044_dp  
    case("Fe")
      mass = 55.845_dp     
    case("Co")
      mass = 58.933194_dp  
    case("Ni")
      mass = 58.6934_dp    
    case("Cu")
      mass = 63.546_dp     
    case("Zn")
      mass = 65.38_dp      
    case("Ga")
      mass = 69.723_dp     
    case("Ge")
      mass = 72.640_dp      
    case("As")
      mass = 74.9216_dp    
    case("Se")
      mass = 78.960_dp     
    case("Br")
      mass = 79.904_dp     
    case("Kr")
      mass = 83.798_dp     
    case("Rb")
      mass = 85.4678_dp    
    case("Sr")
      mass = 87.620_dp      
    case("Y	")
      mass = 88.90584_dp   
    case("Zr")
      mass = 91.224_dp     
    case("Nb")
      mass = 92.90637_dp   
    case("Mo")
      mass = 95.940_dp      
    case("Tc")
      mass = 97.0_dp     
    case("Ru")
      mass = 101.07_dp     
    case("Rh")
      mass = 102.9055_dp   
    case("Pd")
      mass = 106.42_dp     
    case("Ag")
      mass = 107.8682_dp   
    case("Cd")
      mass = 112.414_dp    
    case("In")
      mass = 114.818_dp    
    case("Sn")
      mass = 118.71_dp     
    case("An")
      mass = 121.76_dp     
    case("Te")
      mass = 127.6_dp      
    case("I	")
      mass = 126.90447_dp  
    case("Xe")
      mass = 131.293_dp    
    case("Cs")
      mass = 132.905452_dp
    case("Ba")
      mass = 137.327_dp    
    case("La")
      mass = 138.90547_dp  
    case("Ce")
      mass = 140.116_dp    
    case("Pr")
      mass = 140.90766_dp  
    case("Nd")
      mass = 144.242_dp    
    case("Pm")
      mass = 145.000_dp    
    case("Sa")
      mass = 150.36_dp     
    case("Eu")
      mass = 151.964_dp    
    case("Gd")
      mass = 157.25_dp     
    case("Tb")
      mass = 158.92535_dp  
    case("Dy")
      mass = 162.5_dp      
    case("Ho")
      mass = 164.93033_dp  
    case("Er")
      mass = 167.259_dp    
    case("Tm")
      mass = 168.93422_dp  
    case("Yb")
      mass = 173.045_dp    
    case("Lu")
      mass = 174.9668_dp   
    case("Hf")
      mass = 178.49_dp     
    case("Ta")
      mass = 180.94788_dp  
    case("Tg")
      mass = 183.84_dp     
    case("Re")
      mass = 186.207_dp    
    case("Os")
      mass = 190.23_dp     
    case("Ir")
      mass = 192.217_dp    
    case("Pt")
      mass = 195.084_dp    
    case("Au")
      mass = 196.966569_dp
    case("Hg")
      mass = 200.592_dp
    case("Tl")
      mass = 204.38_dp
    case("Pb")
      mass = 207.2_dp
    case("Bi")
      mass = 208.9804_dp
    case("Po")
      mass = 209.0_dp
    case("At")
      mass = 210.0_dp
    case("Rn")
      mass = 222.0_dp
    case("Fr")
      mass = 223.0_dp
    case("Ra")
      mass = 226.0_dp
    case("Ac")
      mass = 227.0_dp
    case("Th")
      mass = 232.0377_dp
    case("Pa")
      mass = 231.03588_dp
    case("U	")
      mass = 238.02891_dp
    end select 

  end function getAtomicMass

end module dftbp_masses  
