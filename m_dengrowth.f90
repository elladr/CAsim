MODULE m_dengrowth

implicit none
integer, allocatable:: I_v(:,:) 
integer, allocatable:: I_v_G(:,:) 
real, allocatable:: x_G(:,:), y_G(:,:)
real, allocatable:: x_S10(:,:), y_S10(:,:), x_S01(:,:), y_S01(:,:)
real, allocatable:: x_Sn10(:,:), y_Sn10(:,:), x_S0n1(:,:), y_S0n1(:,:)
real, allocatable:: width_cell(:,:), width_cell_old(:,:)
real, allocatable:: V_tip(:,:), V_tip_old(:,:)
real, allocatable:: Teta(:,:)
real, allocatable:: tp(:), tp_init(:), tp_G(:,:)
real, allocatable:: deltaT_old(:), deltaT(:)

END MODULE m_dengrowth