#ifndef _TWALKMISC
#define _TWALKMISC
#include <search.h>
typedef struct node_t
{
    /* Callers expect this to be the first element in the structure - do not
       move!  */
    const void *key;
    struct node_t *left;
    struct node_t *right;
    unsigned int red:1;
} *node;
typedef const struct node_t *const_node;
typedef void (*__action_fn) (const void *, VISIT, int, void * misc);


/* Walk the nodes of a tree.
   ROOT is the root of the tree to be walked, ACTION the function to be
   called at each node.  LEVEL is the level of ROOT in the whole tree.  */
void trecursemisc (const void *vroot, __action_fn action,
        int level, void * misc)
{
    const_node root = (const_node) vroot;

    if (root->left == NULL && root->right == NULL)
        (*action) (root, leaf, level, misc);
    else
    {
        (*action) (root, preorder, level, misc);
        if (root->left != NULL)
            trecursemisc (root->left, action, level + 1 ,misc);
        (*action) (root, postorder, level, misc);
        if (root->right != NULL)
            trecursemisc (root->right, action, level + 1, misc);
        (*action) (root, endorder, level, misc);
    }
}

void twalkmisc (const void *vroot, __action_fn action, void * misc)
{
    const_node root = (const_node) vroot;

    if (root != NULL && action != NULL)
        trecursemisc (root, action, 0, misc);
}

#endif // header guard
