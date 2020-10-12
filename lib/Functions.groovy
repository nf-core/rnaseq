/*
 * This file holds several general utility functions used within the pipeline
 */

class Functions {
    /*
     * This map makes a copy of a map and adds all values specified in a second map
     */
    private static LinkedHashMap map_copy_append(Map template, Map override) {
        def copy = new LinkedHashMap(template)
        copy.putAll(override)
        return copy
    }
}
