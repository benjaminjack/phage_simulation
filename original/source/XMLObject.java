
///////////////////////////////////////////////
//  XMLObject.java
// 
//  @author Jonathan Goler
//  @creation-date 4 August 2001
//   @version 2.0 including GHO mods
//  (c) 2001 Jonathan Goler/MIT Media Lab
//////////////////////////////////////////////


import java.lang.*;
import java.util.*;
import java.io.*;

public class XMLObject {

    public Vector childrenNodes = new Vector();
    public Properties attributes = new Properties();
    public String name = new String();
    public String contents = new String();
    public XMLObject() {};
    /**
     * Main constructor for an XML Object, takes the XML text and creates an XML Object 
     * encapsulating and enabling access to taht data
     *@param xml The XML string that will form the XML Object
     **/
    public XMLObject(String xml) {
	doParse(xml);
    }
    private void doParse(String xml) {
	xml = xml.trim();
	//System.out.println("parse" +xml);
	int i,j;
	int selfCloseEnd = xml.indexOf("/>");
        i = xml.indexOf(">");
	j = xml.indexOf(" ");
	int end = (( i > j) && (j>-1))?j:i;
	//	System.out.println("end="+end);
	//System.out.println(xml.substring(0,end));
	//System.out.println( xml.indexOf(">"));
	name = xml.substring(1,end);
	int closebracket = xml.indexOf(">");
	if(selfCloseEnd == closebracket -1) closebracket--;
	if(end < closebracket) 
	    parseAttributes(xml.substring(end+1, closebracket));
	int childrenStart = xml.indexOf(">")+1;
	int childrenEnd = xml.lastIndexOf("<")-1;
	if(childrenStart<childrenEnd) {
	    if(-1 < name.indexOf("ENCRYPTED")) {
		contents = xml.substring(xml.indexOf(">")+1, xml.lastIndexOf("<"));
		return;
	    }
	    if(childrenStart != xml.indexOf("/>")+2) {
		//	System.out.println("Parsing Children");
		parseChildren(xml.substring(xml.indexOf(">")+1, xml.lastIndexOf("<")));
	    }
	}
    }
    public XMLObject(File file) {
	StringBuffer sb = new StringBuffer();
	try {
	    BufferedReader br = new BufferedReader(new FileReader(file));
	    while(br.ready()) 
		sb.append(br.readLine() + "\n");

	    doParse(sb.toString());
	    
	} catch (Exception e) {
	    System.out.println("Failed to load XML from File");
	    e.printStackTrace();
	}
    }
    public void setContents(String ct) {
	contents = ct;
    }
    public String getContents() {
	return contents;
    }
    /** 
     * Returns the XML children nodes for this object
     **/
    public Vector getChildren() {
	return childrenNodes;
    }
    /**
     * Adds a child node
     *@param child The child node to add 
     **/
    public void addChild(XMLObject child) {
	childrenNodes.addElement(child);
    }
   /**
     * Adds a child node
     *@param child An XML string containing child node to add 
     **/
    public void addChild(String xml) {
	childrenNodes.addElement(new XMLObject(xml));
    }
    /**
     * Returns the attributes of this XML Object
     **/
    public Properties getAttributes() {
	return attributes;
    }
    public void setAttributes(Properties atts) {
	attributes = atts;
    }
    /** 
     * Returns the name of this XML Object
     **/
    public String getName() {
	return name;
    }
    /**
     * Sets the name of this XML Object
     *@param n The name to assign to the object
     **/
    public void setName(String n) {
	name = n;
    }
    /**
     * Sets an Attribute for the XML object
     *@param key the name of the attribute to set
     *@param value the value of the attribute
     **/
    public void setAttribute(String key, String value) {
	attributes.setProperty(key, value);
    }
    public void setProperty(String key, String value) {
	setAttribute(key,value);
    }
    /**
     * Removes an attribute with the key given
     *@param key the name of the attribute to remove
     * 
     **/
    public void removeAttribute(String key) {
	attributes.remove(key);
    }
   
    /**
     * Returns the value of the attribute with the key given
     *@param key the name of the attribute to get
     * 
     **/
    public String getAttribute(String key) {
	return attributes.getProperty(key);
    }
    public String getProperty(String key) {
	return getAttribute(key);
    }
    
    private void parseAttributes(String xml) {
	StringTokenizer tokenizer = new StringTokenizer(xml);
	//	System.out.println("parsing Attributes");
	try{
	while(tokenizer.hasMoreTokens()) {
	    String token = tokenizer.nextToken();
	    String attrName = token.substring(0,token.indexOf("="));
	    String attrValue;
	    if(token.indexOf("\"")>-1) {
		while(token.lastIndexOf("\"")<=token.indexOf("\"")) token = token + " "+tokenizer.nextToken();

			attrValue= token.substring(token.indexOf("\"")+1, 
					   token.lastIndexOf("\""));

	    } else {
		
		attrValue= token.substring(token.indexOf("=")+1);
	    }
	    attributes.setProperty(attrName,attrValue);
	    
	}
	} catch (Exception e) {
	    System.out.println("Exception while parsing attributes"+e+attributes);
	    e.printStackTrace();
	}
    }
    
    private void parseChildren(String xml) {
	try {
	int workIndex = 0;
        while(xml.length()>0) {
	    xml=xml.trim();
	int spc = xml.indexOf(" ");
	int tab = xml.indexOf("\t");
	int ret = xml.indexOf("\n");
	int lr = xml.indexOf("\r");
	int f = xml.indexOf("\f");
	int closing = xml.indexOf(">");
	int selfclosing = xml.indexOf("/>");

	int delimeter = (closing>spc && spc >-1)?spc:closing;

	//	System.out.println("Space="+spc+"t:"+tab+"lf"+ret+ "d"+delimeter+"lr"+lr+"f"+f);
	//\System.out.println(xml.indexOf("<")+1);
	
	if(xml.indexOf("<")< 0) {
		contents = xml;
		return;
	}
	String objName = xml.substring(xml.indexOf("<")+1, delimeter);
	//System.out.println("And here"+objName+">>");	
	String closingName = "</"+objName+">";
	//System.out.println("Even here" + closingName+" index ="+xml.indexOf(closingName));
	if(xml.indexOf(closingName)>0) {
	    //  System.out.println("here");
	    String childXML = xml.substring(0, xml.indexOf(closingName)) +"</"+objName;
	    //System.out.println("but not here");
	  
	    XMLObject temp = new XMLObject(childXML);
	    
	    childrenNodes.add(temp); 
	    
	    xml = xml.substring(childXML.length()+1);
	} else if( selfclosing == closing-1){ 
	 
	    XMLObject temp = new XMLObject(xml.substring(0,closing+1));
	    
	    childrenNodes.add(temp);
	    xml=xml.substring(closing+1);;

	} else xml = "";
	}
	} catch (Exception e) {
	    System.out.println("Exception :"+e);
	    e.printStackTrace();
	    System.out.println("XML="+xml);
	}
    }
/**
 * Returns a string containing the XML text representation of this object
 **/
    public String print() {
	StringBuffer buffer = new StringBuffer();
	buffer.append("<"+name);
        Enumeration props = attributes.propertyNames();
	while(props.hasMoreElements()) {
	    String temp = (String)props.nextElement();
	    temp = temp.replace(' ','_');
	    
	    buffer.append(" "+temp+"=\""+attributes.getProperty(temp)+"\"");
	    
	}

	if(contents.length()==0 && childrenNodes.size()==0){
	    buffer.append("/>\n");
	    return buffer.toString();
	}
	buffer.append(">");
	buffer.append(contents);
	Iterator children = childrenNodes.iterator();
	while(children.hasNext()) { 
	    
	    XMLObject child = (XMLObject)children.next();
	    buffer.append(child.print());
	}
	buffer.append("</"+name+">\n");
	return buffer.toString();
    }
	      
    public XMLObject getChild(String blockHeader) {
	if(name.equals(blockHeader)) return this;
	Iterator i = childrenNodes.iterator();
	while(i.hasNext()) {
	    
	    XMLObject temp = (XMLObject) i.next();
	    XMLObject test = temp.getChild(blockHeader);
	    if(test != null) return test;

	}
	return null;
    }
    /*
     * returns all children of a specified name, so if you have a
     * bunch of nested objects of the same type, say ball, then 
     * the procedure will return all of those children in a vector
     *@param name - String containing the name of the child type
     */
    public Vector getNamedChildren(String name) {
	Vector temp = new Vector();
	for(Iterator i= childrenNodes.iterator(); i.hasNext();) {
	    XMLObject obj = (XMLObject) i.next();
	    if(name.equals(obj.getName())) 
		temp.add(obj);
	}

	return temp;
    }
    public Properties getNVPairs(String attrTag, String nameTag, String valTag) {
	Properties temp = new Properties();
	Vector kiddies = getNamedChildren(attrTag);
	for (Iterator i = kiddies.iterator(); i.hasNext();) {
	    XMLObject child = (XMLObject) i.next();
	    String attrName = child.getAttribute(nameTag);
	    String attrVal = child.getAttribute(valTag);
	    temp.setProperty(attrName, attrVal);

	}

	return temp;
    }

}
    
 
