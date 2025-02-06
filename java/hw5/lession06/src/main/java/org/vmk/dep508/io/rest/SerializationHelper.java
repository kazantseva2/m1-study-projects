package org.vmk.dep508.io.rest;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.log4j.Logger;

import java.io.*;

public class SerializationHelper<T extends Serializable> {

    Class<T> serialazationType;

    public SerializationHelper(Class<T> serialazationType) {
        this.serialazationType = serialazationType;
    }

    private Logger log = Logger.getLogger(getClass());

    ObjectMapper mapper = new ObjectMapper();


    /*
      Необходимо десереализовать объект из файла по указанному пути
     */
    public T loadFromFile(String path) {
        log.info("Attempting to load object from file: " + path);
        try (FileInputStream fileInputStream = new FileInputStream(path);
            ObjectInputStream objectInputStream = new ObjectInputStream(fileInputStream)) {

            T loadedObject = serialazationType.cast(objectInputStream.readObject());
            log.info("Successfully loaded object from file: " + path);
            return loadedObject;

        } catch (IOException | ClassNotFoundException e) {
            log.error("Failed to load object from file: " + path, e);
        }
        return null;
    }

    /*
      Необходимо сохранить сереализованный объект в файл по указанному пути
     */
    public boolean saveToFile(String path, T toSave) {
        log.info("Attempting to save object to file: " + path);
        try (FileOutputStream fileOutputStream = new FileOutputStream(path);
            ObjectOutputStream objectOutputStream = new ObjectOutputStream(fileOutputStream)) {

            objectOutputStream.writeObject(toSave);
            log.info("Successfully saved object to file: " + path);
            return true;

        } catch (IOException e) {
            log.error("Failed to save object to file: " + path, e);
        }
        return false;
    }

    public String convertToJsonString(T toConvert) {
        try {
            String json = mapper.writeValueAsString(toConvert);
            return json;
        } catch (JsonProcessingException e) {
            e.printStackTrace();
        }

        return null;
    }

    public void writeJsonToStream(OutputStream output, T toWrite) {
        try {
            mapper.writeValue(output, toWrite);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public T parseJson(String json) {
        try {
            return mapper.readValue(json, serialazationType);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }
}
