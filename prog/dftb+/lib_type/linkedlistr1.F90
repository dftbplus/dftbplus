#:include 'linkedlist.fypp'

module linkedlistr1
  use accuracy, only : dp
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listRealR1',&
      & ITEM_TYPE='real(dp)',&
      & ITEM_RANK=1,&
      & PADDING='0.0_dp')

end module linkedlistr1
