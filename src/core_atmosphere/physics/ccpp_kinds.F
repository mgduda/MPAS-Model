module ccpp_kinds
!  integer, parameter :: kind_phys = selected_real_kind(6)
   use mpas_kind_types,only: kind_phys => RKIND
   contains
end module ccpp_kinds
