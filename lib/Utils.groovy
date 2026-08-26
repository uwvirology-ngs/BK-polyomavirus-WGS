class Utils {

    /*
     * Adds reference genome information as a secondary meta map in 
     * the form of [ [ meta.id, ... ], [ acc, tag, description ] ]
     * and returns a List of meta + ref_info maps
     */
    public static ArrayList add_ref_info_to_meta(meta, refs_tsv) {
        def new_metadata = new ArrayList<>()
        
        refs_tsv.text.readLines().collect { line ->
            def fields = line.split('\t')
            new_metadata.add(
                [ meta, [ acc: fields[0], tag: fields[1], header: fields[2].toString() ] ]
            )
        }

        return new_metadata
    }

    public static String getAnnotation(String s) {
        String ref = s.substring(s.lastIndexOf('_') + 1)
        return ref
    }

    public static String getGenomicRegion(String s) {
        HashMap map = [
            "PV856408.1": "PV856408.1:1177-2265",
            "PV856409.1": "PV856409.1:1177-2265",
            "PV856429.1": "PV856429.1:1167-2255",
            "PV856444.1": "PV856444.1:1177-2265",
            "PV856452.1": "PV856452.1:1177-2265",
            "PX257737.1": "PX257737.1:1193-2281",
            "PZ129835.1": "PZ129835.1:1470-2558",
            "NC_001538.1": "NC_001538.1:1564-2652"
        ]
        return map[s]
    }

    public static int getCDSLen(String s) {
        HashMap map = [
            "PV856408.1": 1089,
            "PV856409.1": 1089,
            "PV856429.1": 1089,
            "PV856444.1": 1089,
            "PV856452.1": 1089,
            "PX257737.1": 1089,
            "PZ129835.1": 1089,
            "NC_001538.1": 1089
        ]
        return map[s]
    }
}