#:include 'linkedlist.fypp'

module linkedlistmc0
  use accuracy, only : mc
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listCharMc',&
      & ITEM_TYPE='character(mc)',&
      & PADDING='""')

end module linkedlistmc0
