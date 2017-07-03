#:include 'linkedlist.fypp'

module linkedlists0
  use assert
  use xmlf90
  implicit none
  private



  $:define_list(&
      & TYPE_NAME='listString',&
      & ITEM_TYPE='character(len=*)',&
      & NODE_TYPE='type(string)',&
      & PADDING="''", )

end module linkedlists0
