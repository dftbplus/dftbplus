#:include 'linkedlist.fypp'

module linkedlistr0
  use accuracy, only : dp
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listReal',&
      & ITEM_TYPE='real(dp)',&
      & PADDING='0.0_dp')

end module linkedlistr0
