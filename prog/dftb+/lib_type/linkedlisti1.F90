#:include 'linkedlist.fypp'

module linkedlisti1
  use accuracy, only : dp
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listIntR1',&
      & ITEM_TYPE='integer',&
      & ITEM_RANK=1,&
      & PADDING='0')

end module linkedlisti1
