module precision

! Real kinds

integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real

! Integer kinds

integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer

!Complex kinds

integer, parameter :: kc4 = kr4                            ! single precision complex
integer, parameter :: kc8 = kr8                            ! double precision complex

! Text lines
integer, parameter :: mcm_text = 2024
integer, parameter :: rdfpot_text = 2024
integer, parameter :: input_text = 2024
integer, parameter :: filename_text = 64 ! Names of files (RDF, potential, starting, etc.) provided in the main input file

integer, parameter :: atom_name_len = 8 ! Default length for atom/site/atomtype/sitetype name
integer, parameter :: mol_name_len = 12
integer, parameter :: bond_name_len = 20

end module precision
