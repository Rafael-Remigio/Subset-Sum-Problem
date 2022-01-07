// Halll System
typedef struct Heap Heap;
struct Heap{

    // type = 0 (Heap) || type = 1 (MaxHeap)
    char type;
    //
    integer_t *arr;
    // Current Size of the Heap
    integer_t size;
    // Maximum capacity of the heap
    integer_t capacity;
};

integer_t parent(integer_t i){

    // Get the index of the parent
    return (i - 1) / 2;
}

integer_t left_child(integer_t i){

    return (2 * i + 1);
}

integer_t right_child(integer_t i){

    return (2 * i + 2);
}

integer_t get_min(Heap *heap){

    // Only for MinHeaps
    // Return the root node element,
    // since that's the minimum
    return heap->arr[0];
}

integer_t get_max(Heap *heap){

    // Only for MaxHeaps
    // Return the root node element,
    // since that's the maximum
    return heap->arr[0];
}

Heap *init_heap(integer_t capacity){

    Heap *heap = (Heap *)calloc(1ull, sizeof(Heap));
    heap->arr = (integer_t *)calloc(capacity, sizeof(integer_t));
    heap->capacity = capacity;
    heap->size = 0;
    return heap;
}

Heap *insert_heap(Heap *heap, integer_t element, char type){

    if (type == 0){
    
        // Inserts an element to the min heap
        // We first add it to the bottom (last level)
        // of the tree, and keep swapping with it's parent
        // if it is lesser than it. We keep doing that until
        // we reach the root node. So, we will have inserted the
        // element in it's proper position to preserve the min heap property
        if (heap->size == heap->capacity)
        {
            fprintf(stderr, "Cannot insert %llu. Heap is already full!\n", element);
            return heap;
        }
        // We can add it. Increase the size and add it to the end
        heap->size++;
        heap->arr[heap->size - 1ull] = element;

        // Keep swapping until we reach the root
        integer_t curr = heap->size - 1ull;
        // As long as you aren't in the root node, and while the
        // parent of the last element is greater than it
        while (curr > 0 && heap->arr[parent(curr)] > heap->arr[curr])
        {
            // Swap
            integer_t temp = heap->arr[parent(curr)];
            heap->arr[parent(curr)] = heap->arr[curr];
            heap->arr[curr] = temp;
            // Update the current index of element
            curr = parent(curr);
        }
        return heap;
    }
    else if (type == 1){
    
        // Inserts an element to the max heap
        // We first add it to the bottom (last level)
        // of the tree, and keep swapping with it's parent
        // if it is grater than it. We keep doing that until
        // we reach the root node. So, we will have inserted the
        // element in it's proper position to preserve the max heap property
        if (heap->size == heap->capacity)
        {
            fprintf(stderr, "Cannot insert %llu. Heap is already full!\n", element);
            return heap;
        }
        // We can add it. Increase the size and add it to the end
        heap->size++;
        heap->arr[heap->size - 1ull] = element;

        // Keep swapping until we reach the root
        integer_t curr = heap->size - 1ull;
        // As long as you aren't in the root node, and while the
        // parent of the last element is lesser than it
        while (curr > 0 && heap->arr[parent(curr)] < heap->arr[curr])
        {
            // Swap
            integer_t temp = heap->arr[parent(curr)];
            heap->arr[parent(curr)] = heap->arr[curr];
            heap->arr[curr] = temp;
            // Update the current index of element
            curr = parent(curr);
        }
        return heap;
    }
    return heap;
}

Heap *heapify(Heap *heap, integer_t index, char type){

    switch (type)
    {
    case 0:
        // Rearranges the heap as to maintain
        // the min-heap property
        if (heap->size <= 1ull)
            return heap;

        integer_t left = left_child(index);
        integer_t right = right_child(index);

        // Variable to get the smallest element of the subtree
        // of an element an index
        integer_t smallest = index;

        // If the left child is smaller than this element, it is
        // the smallest
        if (left < heap->size && heap->arr[left] < heap->arr[index])
            smallest = left;

        // Similarly for the right, but we are updating the smallest element
        // so that it will definitely give the least element of the subtree
        if (right < heap->size && heap->arr[right] < heap->arr[smallest])
            smallest = right;

        // Now if the current element is not the smallest,
        // swap with the current element. The min heap property
        // is now satisfied for this subtree. We now need to
        // recursively keep doing this until we reach the root node,
        // the point at which there will be no change!
        if (smallest != index)
        {
            integer_t temp = heap->arr[index];
            heap->arr[index] = heap->arr[smallest];
            heap->arr[smallest] = temp;
            heap = heapify(heap, smallest, 0);
        }

        return heap;
        break;
    case 1:
        // Rearranges the heap as to maintain
        // the max-heap property
        if (heap->size <= 1ull)
            return heap;

        left = left_child(index);
        right = right_child(index);

        // Variable to get the greatest element of the subtree
        // of an element an index
        integer_t greatest = index;

        // If the left child is greatest than this element, it is
        // the greatest
        if (left < heap->size && heap->arr[left] > heap->arr[index])
            greatest = left;

        // Similarly for the right, but we are updating the greatest element
        // so that it will definitely give the least element of the subtree
        if (right < heap->size && heap->arr[right] > heap->arr[greatest])
            greatest = right;

        // Now if the current element is not the greatest,
        // swap with the current element. The max heap property
        // is now satisfied for this subtree. We now need to
        // recursively keep doing this until we reach the root node,
        // the point at which there will be no change!
        if (greatest != index)
        {
            integer_t temp = heap->arr[index];
            heap->arr[index] = heap->arr[greatest];
            heap->arr[greatest] = temp;
            heap = heapify(heap, greatest, 1);
        }

        return heap;
        break;
    }
    return heap;
}

Heap *delete_minimum(Heap *heap){

    // Deletes the minimum element, at the root
    if (!heap || heap->size == 0)
        return heap;

    integer_t size = heap->size;
    integer_t last_element = heap->arr[size - 1];

    // Update root value with the last element
    heap->arr[0] = last_element;

    // Now remove the last element, by decreasing the size
    heap->size--;
    size--;

    // We need to call heapify(), to maintain the min-heap
    // property
    heap = heapify(heap, 0, 0);
    return heap;
}

Heap *delete_maximum(Heap *heap){

    // Deletes the maximum element, at the root
    if (!heap || heap->size == 0)
        return heap;

    integer_t size = heap->size;
    integer_t last_element = heap->arr[size - 1];

    // Update root value with the last element
    heap->arr[0] = last_element;

    // Now remove the last element, by decreasing the size
    heap->size--;
    size--;

    // We need to call heapify(), to maintain the max-heap
    // property
    heap = heapify(heap, 0, 1);
    return heap;
}

Heap *delete_element(Heap *heap, integer_t index, char type){

    switch (type)
    {
    case 0:
        // Deletes an element, indexed by index
        // Ensure that it's lesser than the current root
        heap->arr[index] = get_min(heap) - 1;

        // Now keep swapping, until we update the tree
        integer_t curr = index;
        while (curr > 0 && heap->arr[parent(curr)] > heap->arr[curr])
        {
            integer_t temp = heap->arr[parent(curr)];
            heap->arr[parent(curr)] = heap->arr[curr];
            heap->arr[curr] = temp;
            curr = parent(curr);
        }

        // Now simply delete the minimum element
        heap = delete_minimum(heap);
        return heap;
        break;
    case 1:
        // Deletes an element, indexed by index
        // Ensure that it's lesser than the current root
        heap->arr[index] = get_min(heap) - 1;

        // Now keep swapping, until we update the tree
        curr = index;
        while (curr > 0 && heap->arr[parent(curr)] < heap->arr[curr])
        {
            integer_t temp = heap->arr[parent(curr)];
            heap->arr[parent(curr)] = heap->arr[curr];
            heap->arr[curr] = temp;
            curr = parent(curr);
        }

        // Now simply delete the maximum element
        heap = delete_maximum(heap);
        return heap;
        break;
        break;
    }
    return heap;
}