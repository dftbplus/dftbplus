#:include 'linkedlist.fypp'

module linkedlistr2
  use accuracy, only : dp
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listRealR2',&
      & ITEM_TYPE='real(dp)',&
      & ITEM_RANK=2,&
      & PADDING='0.0_dp')

end module linkedlistr2
