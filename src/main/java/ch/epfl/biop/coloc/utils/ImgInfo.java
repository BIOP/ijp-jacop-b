/*-
 * #%L
 * An update to the JACoP Plugin that helps in the management of ROIs, Z sections and helps generate cleaner reports
 * %%
 * Copyright (C) 2009 - 2026 Susanne Bolte, Fabrice P. Cordelières, Olivier Burri
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
/*
 * ImgInfo.java
 *
 * Created on 14 janvier 2008, 21:11
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package ch.epfl.biop.coloc.utils;


/**
 *
 * @author Fabrice Cordeli�res
 */
public class ImgInfo {
    public String title;
    public int min;
    public int max;
    public int thr;
    
    /** Creates a new instance of ImgInfo */
    public ImgInfo() {
        this.title="[No image]";
        this.min=0;
        this.max=0;
        this.thr=0;
    }
    /** Creates a new instance of ImgInfo */
    public ImgInfo(String title, int min, int max, int thr) {
        this.title=title;
        this.min=min;
        this.max=max;
        this.thr=thr;
     }
    
}
