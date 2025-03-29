def element_to_dict(element):
    """
    Convert an ElementTree Element object to a dictionary.
    """
    result = {}
    result.update(element.attrib)

    # Process the element's children recursively
    for child in element:
        # Recursively convert the child element to a dictionary
        child_dict = element_to_dict(child)
        # If the child tag already exists in the dictionary, convert it to a list
        if child.tag in result:
            if not isinstance(result[child.tag], list):
                result[child.tag] = [result[child.tag]]
            result[child.tag].append(child_dict)
        else:
            result[child.tag] = child_dict

    # If the element has text, add it to the dictionary after stripping whitespace
    if element.text:
        result[element.tag] = element.text.strip()

    def remove_namespace(tag):
        return tag.split('}', 1)[-1] if '}' in tag else tag
    
    def clean_dict(d):
        if isinstance(d, dict):
            return {remove_namespace(k): clean_dict(v) for k, v in d.items()}
        elif isinstance(d, list):
            return [clean_dict(i) for i in d]
        else:
            return d
    
    result = clean_dict(result)
    return result

def rel_elsize(elsize):
    """
    Calculate the relative element sizes for an array, smallest is set to 1
    """
    min_size=min(elsize)
    elsize_scaled= [round(x/min_size, 2) for x in elsize]
    return(elsize_scaled)