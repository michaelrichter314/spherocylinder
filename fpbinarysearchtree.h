#pragma once

// item in a binary tree. data is the payload, rest are pointers to other elements
template<typename type>
struct BinaryTreeItem
{
	type data;
	BinaryTreeItem* rightChild = NULL;
	BinaryTreeItem* leftChild = NULL;
	BinaryTreeItem* parent = NULL;
};

// binary search tree with internal storage of items
// "type" must be sortable 

template <typename type>
class BinarySearchTree
{
	public:

		// destructor. necessary to free memory upon destruction
		~BinarySearchTree();

		// inser an item
		void insert(type);

		// remove an item from tree, free memory associated with item
		void remove(type);

		// remove all items
		void clear();

		// start iterating through items (sorted)
		void startTraversal();

		// get next value in iteration. undefined behavior if no more elements or never called startTraversal
		type getNextElement();

		// getNextElement possible?
		bool gotAnotherElement() const;

		// returns lowest element
		type getLowestElement();

		// returns highest element
		type getHighestElement();

		//void printElements();
		//void printElements(BinaryTreeItem<type>*);

	private:
		// root node
		BinaryTreeItem<type>* root = NULL;

		// current position in iteration through elements
		BinaryTreeItem<type>* currentItem = NULL;

		BinaryTreeItem<type>* getLowestItem();
		BinaryTreeItem<type>* getHighestItem();

		// remove one item, fix sourrounding pointers
		void remove(BinaryTreeItem<type>*);

		// remove item, do not fix pointers
		void removeNoUpdate(BinaryTreeItem<type>*);

		// get item by value
		BinaryTreeItem<type>* getItem(type);

		// returns pointer to left most grand-grand-... child
		BinaryTreeItem<type>* getLeftmostChild(BinaryTreeItem<type>*);

		// removes node + all subnodes, does not fix pointers
		void removeSubtreeNoUpdate(BinaryTreeItem<type>*);

};

template<typename type>
BinarySearchTree<type>::~BinarySearchTree()
{
	clear();
}

/*template<typename type>
void BinarySearchTree<type>::printElements()
{
	printElements(root);
}
template<typename type>
void BinarySearchTree<type>::printElements(BinaryTreeItem<type>* item)
{
	if (item == NULL) {return;}
	printElements(item -> leftChild);
	std::cout << item -> data << std::endl;
	printElements(item -> rightChild);
}*/

template<typename type>
void BinarySearchTree<type>::insert(type data)
{
	BinaryTreeItem<type>* newItem = new(BinaryTreeItem<type>);
	newItem -> data = data;

	if (root == NULL)
	{
		root = newItem;
	}
	else
	{
		BinaryTreeItem<type>* parent = root;
		while (true) // will terminate unless data corrupted
		{
			if (data <= (parent->data))
			{
				if (parent -> leftChild != NULL)
				{
					parent = parent -> leftChild;
				}
				else
				{
					parent -> leftChild = newItem;
					newItem -> parent = parent;
					break;
				}
			}
			else
			{
				if (parent -> rightChild != NULL)
				{
					parent = parent -> rightChild;
				}
				else
				{
					parent -> rightChild = newItem;
					newItem -> parent = parent;
					break;
				}
			}
		}
	}
}


template<typename type>
void BinarySearchTree<type>::remove(type data)
{
	auto item = getItem(data);
	if (item != NULL)
	{
		remove(item);
	}
}

template<typename type>
void BinarySearchTree<type>::clear()
{
	removeSubtreeNoUpdate(root);
}

template<typename type>
void BinarySearchTree<type>::startTraversal()
{
	currentItem = getLeftmostChild(root);
}

template<typename type>
bool BinarySearchTree<type>::gotAnotherElement() const
{
	return (currentItem != NULL);
}

template<typename type>
type BinarySearchTree<type>::getNextElement()
{
	// currentItem is already set. save return value
	type toReturn = currentItem -> data;

	// update currentItem
	if (currentItem -> rightChild != NULL)
	{
		currentItem = getLeftmostChild(currentItem -> rightChild);
	}
	else
	{
		// back up one step
		bool haveToBackUpAgain = true;
		while (haveToBackUpAgain == true && currentItem -> parent != NULL)
		{
			// check if we are going to upper left or right
			if (currentItem -> parent -> rightChild == currentItem)
			{
				haveToBackUpAgain = true;
			}
			else
			{
				haveToBackUpAgain = false;
			}
			currentItem = currentItem -> parent;

		}
			if (haveToBackUpAgain) {currentItem = NULL;}
	}

	return toReturn;
}





template<typename type>
BinaryTreeItem<type>* BinarySearchTree<type>::getItem(type data)
{
	auto item = root;

	while (item != NULL && item->data != data)
	{
		if (data <= (item->data))
		{
			item = item -> leftChild;
		}
		else
		{
			item = item -> rightChild;
		}
	}

	return item;
}

template<typename type>
void BinarySearchTree<type>::removeSubtreeNoUpdate(BinaryTreeItem<type>* item)
{
	if (item == NULL) {return;}

	removeSubtreeNoUpdate(item -> leftChild);
	removeSubtreeNoUpdate(item -> rightChild);
	removeNoUpdate(item);
}

template<typename type>
void BinarySearchTree<type>::removeNoUpdate(BinaryTreeItem<type>* item)
{
	if (item == NULL) {return;}

	delete(item);
}

template<typename type>
void BinarySearchTree<type>::remove(BinaryTreeItem<type>* item)
{
	if (item == NULL) {return;}

	int children = 0;
	if (item -> leftChild != NULL) {children++;}
	if (item -> rightChild != NULL) {children++;}

	if (children == 0)
	{		
		// no children. just fix parent
		if (item -> parent != NULL)
		{
			if (item -> parent -> leftChild == item)
			{
				item -> parent -> leftChild = NULL;
			}
			else
			{
				item -> parent -> rightChild = NULL;
			}
		}
		else
		{
			root = NULL;
		}
		delete(item);
	}
	else if (children == 1)
	{
		// one child. stitch child with parent
		if (item -> leftChild != NULL)
		{
			item -> leftChild -> parent = item -> parent;
			if (item -> parent != NULL)
			{
				if (item -> parent -> leftChild == item)
				{
					item -> parent -> leftChild = item -> leftChild;
				}
				else
				{
					item -> parent -> rightChild = item -> leftChild;
				}
			}
			else
			{
				root = item -> leftChild;
			}
		}
		else
		{
			item -> rightChild -> parent = item -> parent;
			if (item -> parent != NULL)
			{
				if (item -> parent -> leftChild == item)
				{
					item -> parent -> leftChild = item -> rightChild;
				}
				else
				{
					item -> parent -> rightChild = item -> rightChild;
				}
			}
			else
			{
				root = item -> rightChild;
			}

		}
		delete(item);
	}
	else
	{
		auto leftEnd = getLeftmostChild(item->rightChild);
		item -> data = leftEnd -> data;

		remove(leftEnd);
	}	

}
template <typename type>
BinaryTreeItem<type>* BinarySearchTree<type>::getLeftmostChild(BinaryTreeItem<type>* item)
{
	if (item == NULL) {return NULL;}

	while (item -> leftChild != NULL) {item = item -> leftChild;}

	return item;
}

