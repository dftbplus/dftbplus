#:include 'linkedlist.fypp'

module linkedlistlc0
  use accuracy, only : lc
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listCharLc',&
      & ITEM_TYPE='character(lc)',&
      & PADDING='""')

end module linkedlistlc0
